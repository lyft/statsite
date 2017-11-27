#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <assert.h>
#include <pthread.h>
#include <signal.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <math.h>

#include "likely.h"
#include "metrics.h"
#include "sink.h"
#include "streaming.h"
#include "conn_handler.h"
#include "config.h"

/* Static method declarations */
static int handle_ascii_client_connect(statsite_conn_handler *handle);
static int buffer_after_terminator(char *buf, int buf_len, char terminator, char **after_term, int *after_len);

/**
 * This is the current metrics object we are using
 */
static metrics *GLOBAL_METRICS;
static statsite_config *GLOBAL_CONFIG;

/**
 * Invoked to initialize the conn handler layer.
 */
void init_conn_handler(statsite_config *config) {
    // Make the initial metrics object
    metrics *m = malloc(sizeof(metrics));
    int res = init_metrics(config->timer_eps, config->tdigest_compression, config->quantiles,
            config->num_quantiles, config->histograms, config->set_precision, m);
    assert(res == 0);
    GLOBAL_METRICS = m;

    // Store the config
    GLOBAL_CONFIG = config;
}

/**
 * A struct passed to the flush thread which contains an instance of
 * metrics and any currently configured sinks.
 */
struct flush_op {
    metrics* m;
    sink* sinks;
};

/**
 * This is the thread that is invoked to handle flushing metrics
 */
static void* flush_thread(void *arg) {
    // Cast the args
    struct flush_op* ops = arg;
    metrics *m = ops->m;
    sink* sinks = ops->sinks;

    // Get the current time
    struct timeval tv;
    gettimeofday(&tv, NULL);

    for (sink* s = sinks; s != NULL; s = s->next) {
        int res = s->command(s, m, &tv);
        if (res != 0) {
            syslog(LOG_WARNING, "Streaming command %s exited with status %d", s->sink_config->name, res);
        }
    }

    // Cleanup
    destroy_metrics(m);
    free(m);
    free(ops);
    return NULL;
}

/**
 * Invoked to when we've reached the flush interval timeout
 */
void flush_interval_trigger(sink* sinks) {
    // Make a new metrics object
    metrics *m = malloc(sizeof(metrics));
    init_metrics(GLOBAL_CONFIG->timer_eps, GLOBAL_CONFIG->tdigest_compression, GLOBAL_CONFIG->quantiles,
            GLOBAL_CONFIG->num_quantiles, GLOBAL_CONFIG->histograms,
            GLOBAL_CONFIG->set_precision, m);

    // Swap with the new one
    struct flush_op* ops = calloc(1, sizeof(struct flush_op));
    ops->m = GLOBAL_METRICS;
    ops->sinks = sinks;

    GLOBAL_METRICS = m;

    // Start a flush thread
    pthread_t thread;
    sigset_t oldset;
    sigset_t newset;
    sigfillset(&newset);
    pthread_sigmask(SIG_BLOCK, &newset, &oldset);
    int err = pthread_create(&thread, NULL, flush_thread, ops);
    pthread_sigmask(SIG_SETMASK, &oldset, NULL);

    if (err == 0) {
        pthread_detach(thread);
        return;
    }

    syslog(LOG_WARNING, "Failed to spawn flush thread: %s", strerror(err));
    GLOBAL_METRICS = ops->m;
    destroy_metrics(m);
    free(m);
}

/**
 * Called when statsite is terminating to flush the
 * final set of metrics
 */
void final_flush(sink* sinks) {
    // Get the last set of metrics
    metrics *old = GLOBAL_METRICS;
    GLOBAL_METRICS = NULL;

    /* We heap allocate this in order to allow it to be freed by the function */
    struct flush_op* ops = calloc(1, sizeof(struct flush_op));
    ops->m = old;
    ops->sinks = sinks;

    flush_thread(ops);

    for (sink* sink = sinks; sink != NULL; sink = sink->next) {
        if (sink->close)
            sink->close(sink);
    }
}


/**
 * Invoked by the networking layer when there is new
 * data to be handled. The connection handler should
 * consume all the input possible, and generate responses
 * to all requests.
 * @arg handle The connection related information
 * @return 0 on success.
 */
int handle_client_connect(statsite_conn_handler *handle) {
    // Try to read the magic character, bail if no data
    unsigned char magic;
    if (unlikely(peek_client_byte(handle->conn, &magic) == -1)) return 0;

    return handle_ascii_client_connect(handle);
}

/**
 * Invoked to handle ASCII commands. This is the default
 * mode for statsite, to be backwards compatible with statsd
 * @arg handle The connection related information
 * @return 0 on success.
 */
static int handle_ascii_client_connect(statsite_conn_handler *handle) {
    // Look for the next command line
    char *buf, *key, *val_str, *type_str, *sample_str, *endptr;
    metric_type type;
    int buf_len, should_free, status, i, after_len;
    double val;
    double sample_rate = 1.0;

    while (1) {
        sample_rate = 1.0;

        status = extract_to_terminator(handle->conn, '\n', &buf, &buf_len, &should_free);
        if (status == -1) return 0; // Return if no command is available

        // Check for a valid metric
        // Scan for the colon
        status = buffer_after_terminator(buf, buf_len, ':', &val_str, &after_len);
        if (likely(!status)) status |= buffer_after_terminator(val_str, after_len, '|', &type_str, &after_len);
        if (unlikely(status)) {
            syslog(LOG_WARNING, "Failed parse metric! Input: %s", buf);
            goto ERR_RET;
        }

        // Convert the type
        switch (*type_str) {
            case 'c':
                type = COUNTER;
                break;
            case 'h':
            case 'm':
                type = TIMER;
                break;
            case 'k':
                type = KEY_VAL;
                break;
            case 'g':
                type = GAUGE;

                // Check if this is a delta update
                switch (*val_str) {
                    case '+':
                    case '-':
                        type = GAUGE_DELTA;
                }
                break;
            case 's':
                type = SET;
                break;
            default:
                type = UNKNOWN;
                syslog(LOG_WARNING, "Received unknown metric type! Input: %c", *type_str);
                goto ERR_RET;
        }

        // Increment the number of inputs received
        if (GLOBAL_CONFIG->input_counter)
            metrics_add_sample(GLOBAL_METRICS, COUNTER, GLOBAL_CONFIG->input_counter, 1, sample_rate);

        // Fast track the set-updates
        if (type == SET) {
            metrics_set_update(GLOBAL_METRICS, buf, val_str);
            goto END_LOOP;
        }

        // Convert the value to a double
        val = strtod(val_str, &endptr);
        if (unlikely(endptr == val_str)) {
            syslog(LOG_WARNING, "Failed value conversion! Input: %s", val_str);
            goto ERR_RET;
        }

        // Handle counter sampling if applicable
        if ((type == COUNTER || type == TIMER) && !buffer_after_terminator(type_str, after_len, '@', &sample_str, &after_len)) {
            double unchecked_rate = strtod(sample_str, &endptr);
            if (unlikely(endptr == sample_str)) {
                syslog(LOG_WARNING, "Failed sample rate conversion! Input: %s", sample_str);
                goto ERR_RET;
            }
            if (likely(unchecked_rate > 0 && unchecked_rate <= 1)) {
                sample_rate = unchecked_rate;
                // Magnify the value
                if (type == COUNTER) {
                    val = val * (1.0 / sample_rate);
                }
            }
        }

        // Store the sample
        metrics_add_sample(GLOBAL_METRICS, type, buf, val, sample_rate);

END_LOOP:
        // Make sure to free the command buffer if we need to
        if (should_free) free(buf);
    }

    return 0;
ERR_RET:
    if (should_free) free(buf);
    return -1;
}

/**
 * Scans the input buffer of a given length up to a terminator.
 * Then sets the start of the buffer after the terminator including
 * the length of the after buffer.
 * @arg buf The input buffer
 * @arg buf_len The length of the input buffer
 * @arg terminator The terminator to scan to. Replaced with the null terminator.
 * @arg after_term Output. Set to the byte after the terminator.
 * @arg after_len Output. Set to the length of the output buffer.
 * @return 0 if terminator found. -1 otherwise.
 */
static int buffer_after_terminator(char *buf, int buf_len, char terminator, char **after_term, int *after_len) {
    // Scan for a space
    char *term_addr = memchr(buf, terminator, buf_len);
    if (!term_addr) {
        *after_term = NULL;
        return -1;
    }

    // Convert the space to a null-seperator
    *term_addr = '\0';

    // Provide the arg buffer, and arg_len
    *after_term = term_addr+1;
    *after_len = buf_len - (term_addr - buf + 1);
    return 0;
}
