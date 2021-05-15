#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <jansson.h>
#include <curl/curl.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <zstd.h>

#include "lifoq.h"
#include "metrics.h"
#include "sink.h"
#include "strbuf.h"
#include "utils.h"
#include "rand.h"
#include "elide.h"

const int MAX_BODY_OBJECTS = 10000;
const int DEFAULT_WORKERS = 2;
const useconds_t FAILURE_WAIT = 5000000; /* 5 seconds */

const char* DEFAULT_CIPHERS_NSS = "ecdhe_ecdsa_aes_128_gcm_sha_256,ecdhe_rsa_aes_256_sha,rsa_aes_128_gcm_sha_256,rsa_aes_256_sha,rsa_aes_128_sha";
const char* DEFAULT_CIPHERS_OPENSSL = "EECDH+AESGCM:EDH+AESGCM:AES256+EECDH:AES256+EDH:DHE-RSA-AES256-SHA:DHE-RSA-AES128-SHA";

const char* USERAGENT = "statsite-http/0";
const char* OAUTH2_GRANT = "grant_type=client_credentials";

struct http_queue_entry {
    time_t not_before_backoff;
    void* data;
};

struct http_sink {
    sink sink;
    lifoq* queue;
    pthread_t* worker;
    pthread_mutex_t sink_mutex;
    char* oauth_bearer;
    elide_t *elide;
    int elide_skip;
};

/*
 * Data from the metrics_iter callback */
struct cb_info {
    json_t** jobjects;
    int jobjects_count;

    elide_t *elide;
    const statsite_config* config;
    const sink_config_http *httpconfig;
    struct timeval now;
};

/**
 * Helper function to return a random double
 */
static double _get_random(void) {
    double random_delay;
    long int rand_seed;
    if (rand_gather((char*)&rand_seed, sizeof(long int)) == -1) {
        syslog(LOG_NOTICE, "Falling back on system time to seed HTTP rand");
        rand_seed = time(0);
    }
    struct drand48_data randbuf;
    srand48_r(rand_seed, &randbuf);
    drand48_r(&randbuf, &random_delay);
    return random_delay;
}

/*
 * Check if the sink's elision output buffer is not initialized, or
 * if it requires re-initialization to purge old metrics.
 */
static void sink_elide_refresh(struct http_sink* s) {
    struct timeval now;
    const sink_config_http* httpconfig = (sink_config_http*)s->sink.sink_config;

    if (httpconfig->elide_interval == 0) {
        return;
    }

    gettimeofday(&now, NULL);
    now.tv_sec -= 60*15;
    if (s->elide == NULL) {
        elide_init(&s->elide, s->elide_skip % httpconfig->elide_interval);
    } else {
        int removed = elide_gc(s->elide, now);
        syslog(LOG_DEBUG, "HTTP: GC elide removed %d entries", removed);
    }
}

/**
 * Check and report if this metric is elided or not.
 * Returns: 0 if not, 1 if elided and should be skipped
 */
static int check_elide(struct cb_info* info, char* full_name, double value) {
    if (info->httpconfig->elide_interval == 0) {
        return 0;
    }

    if (value == 0) {
        int res = elide_mark(info->elide, full_name, info->now);
        if (res % info->httpconfig->elide_interval != info->elide->skip) {
            return 1;
        }
    } else {
        elide_unmark(info->elide, full_name, info->now);
    }
    return 0;
}

/*
 * Callback handling add metrics. If a destination object is too full
 * the callback will work to split the destination object and generate
 * a new one.
 */
static int add_metrics(void* data,
                       metric_type type,
                       char* name,
                       void* value) {
    struct cb_info* info = (struct cb_info*)data;

    if (json_object_size(info->jobjects[0]) > MAX_BODY_OBJECTS) {
        /* build an array size + 1 and copy in references */
        json_t** jobjects = malloc(sizeof(json_t*) * (info->jobjects_count + 1));
        memcpy(&jobjects[1], &info->jobjects[0], sizeof(json_t*) * info->jobjects_count);
        free(info->jobjects);

        info->jobjects = jobjects;
        info->jobjects[0] = json_object();
        info->jobjects_count++;
    }
    json_t* obj = info->jobjects[0];
    const statsite_config* config = info->config;

    /*
     * Scary looking macro to apply suffixes to a string, and then
     * insert them into a json object with a given value. Needs "suffixed"
     * in scope, with a "base_len" telling it how mmany chars the string is
     */
#define SUFFIX_ADD(suf, val)                                            \
    do {                                                                \
        suffixed[base_len - 1] = '\0';                                  \
        strcat(suffixed, suf);                                          \
        json_object_set_new(obj, suffixed, val);                        \
    } while(0)

    char* prefix = config->prefixes_final[type];
    /* Using C99 stack allocation, don't panic */
    int base_len = strlen(name) + strlen(prefix) + 1;
    char full_name[base_len];
    strcpy(full_name, prefix);
    strcat(full_name, name);
    switch (type) {
    case GAUGE_DIRECT:
    {
        double gv = gauge_direct_value(value);
        if (check_elide(info, full_name, gv) == 1)
            break;
        json_object_set_new(obj, full_name, json_real(gauge_direct_value(value)));
        break;
    }
    case GAUGE:
    {
        const int suffix_space = 8;
        char suffixed[base_len + suffix_space];
        double sum = gauge_sum(value);
        strcpy(suffixed, full_name);
        if (check_elide(info, full_name, sum) == 1)
            break;
        json_object_set_new(obj, full_name, json_real(gauge_value(value)));
        SUFFIX_ADD(".sum", json_real(sum));
        SUFFIX_ADD(".mean", json_real(gauge_mean(value)));
        SUFFIX_ADD(".min", json_real(gauge_min(value)));
        SUFFIX_ADD(".max", json_real(gauge_max(value)));
        break;
    }
    case COUNTER:
    {
        if (config->extended_counters) {
            /* We allow up to 8 characters for a suffix, based on the static strings below */
            const int suffix_space = 8;
            char suffixed[base_len + suffix_space];
            strcpy(suffixed, full_name);
            double sum = counter_sum(value);
            if (check_elide(info, full_name, sum) == 1)
                break;
            SUFFIX_ADD(".count", json_integer(counter_count(value)));
            SUFFIX_ADD(".sum", json_real(sum));
            SUFFIX_ADD(".lower", json_real(counter_min(value)));
            SUFFIX_ADD(".upper", json_real(counter_max(value)));
            SUFFIX_ADD(".rate", json_real(counter_sum(value) / config->flush_interval));
        } else {
            json_object_set_new(obj, full_name, json_real(counter_sum(value)));
        }
        break;
    }
    case SET:
        json_object_set_new(obj, full_name, json_integer(set_size(value)));
        break;
    case TIMER:
    {
        timer_hist *t = (timer_hist*)value;

        double mean = timer_mean(&t->tm);
        if (check_elide(info, full_name, mean) == 1)
            break;

        /* We allow up to 40 characters for the metric name suffix. */
        const int suffix_space = 100;
        char suffixed[base_len + suffix_space];
        strcpy(suffixed, full_name);
        SUFFIX_ADD(".mean", json_real(mean));
        SUFFIX_ADD(".lower", json_real(timer_min(&t->tm)));
        SUFFIX_ADD(".upper", json_real(timer_max(&t->tm)));
        SUFFIX_ADD(".count", json_integer(timer_count(&t->tm)));
        for (int i = 0; i < config->num_quantiles; i++) {
            char ptile[suffix_space];
            int percentile;
            double quantile = config->quantiles[i];
            /**
             * config.c already does sanity checks
             * on the quantiles input, dont need to
             * worry about it here.
             */
            to_percentile(quantile, &percentile);
            snprintf(ptile, suffix_space, ".p%d", percentile);
            ptile[suffix_space-1] = '\0';
            SUFFIX_ADD(ptile, json_real(timer_query(&t->tm, quantile)));
        }
        SUFFIX_ADD(".rate", json_real(timer_sum(&t->tm) / config->flush_interval));

        /* Manual histogram bins */
        if (t->conf) {
            char ptile[suffix_space];
            snprintf(ptile, suffix_space, ".bin_<%0.2f", t->conf->min_val);
            ptile[suffix_space-1] = '\0';
            SUFFIX_ADD(ptile, json_integer(t->counts[0]));
            for (int i = 0; i < t->conf->num_bins - 2; i++) {
                sprintf(ptile, ".bin_%0.2f", t->conf->min_val+(t->conf->bin_width*i));
                SUFFIX_ADD(ptile, json_integer(t->counts[i+1]));
            }
            sprintf(ptile, ".bin_>%0.2f", t->conf->max_val);
            SUFFIX_ADD(ptile, json_integer(t->counts[t->conf->num_bins - 1]));
        }
        break;
    }
    default:
        syslog(LOG_ERR, "Unknown metric type: %d", type);
        break;
    }

    return 0;
}

static int json_cb(const char* buf, size_t size, void* d) {
    strbuf* sbuf = (strbuf*)d;
    strbuf_cat(sbuf, buf, size);
    return 0;
}

static int serialize_jobject(struct http_sink* sink,
                            json_t* jobject,
                            struct timeval* tv,
                            time_t not_before) {
    const sink_config_http* httpconfig = (const sink_config_http*)sink->sink.sink_config;

    size_t obj_size = json_object_size(jobject);
    /* Cowardly refuse to make an empty json list */
    if (obj_size == 0) {
        json_decref(jobject);
        return 0;
    }
    strbuf* json_buf;
    strbuf_new(&json_buf, 0);
    json_dump_callback(jobject, json_cb, (void*)json_buf, 0);

    int json_len = 0;
    char* json_data = strbuf_get(json_buf, &json_len);
    json_decref(jobject);

    /* Many APIs reject empty metrics lists - in this case, we simply early abort */
    if (json_len == 2) {
        strbuf_free(json_buf, true);
        return 0;
    }

    /* Start working on the post buffer contents */
    strbuf* post_buf;
    strbuf_new(&post_buf, 128);

    CURL* curl = curl_easy_init();

    /* Encode the json document as a parameter */
    char *escaped_json_data = curl_easy_escape(curl, json_data, json_len);
    strbuf_catsprintf(post_buf, "%s=%s", httpconfig->metrics_name, escaped_json_data);
    strbuf_free(json_buf, true);
    curl_free(escaped_json_data);

    /* Encode the time stamp */
    struct tm tm;
    localtime_r(&tv->tv_sec, &tm);
    char time_buf[200];
    strftime(time_buf, 200, httpconfig->timestamp_format, &tm);
    char* encoded_time = curl_easy_escape(curl, time_buf, 0);
    strbuf_catsprintf(post_buf, "&%s=%s", httpconfig->timestamp_name, encoded_time);
    curl_free(encoded_time);

    /* Encode all the free-form parameters from configuration */
    for (kv_config* kv = httpconfig->params; kv != NULL; kv = kv->next) {
        char* encoded = curl_easy_escape(curl, kv->v, 0);
        strbuf_catsprintf(post_buf, "&%s=%s", kv->k, encoded);
        curl_free(encoded);
    }

    curl_easy_cleanup(curl);

    int post_len = 0;
    char* post_data = strbuf_get(post_buf, &post_len);
    strbuf_free(post_buf, false);

    // Compres the body before enqueueing
    size_t bound = ZSTD_compressBound(post_len);
    void* compressed_post_data = malloc(bound * sizeof(char));
    size_t compressed_size = ZSTD_compress(compressed_post_data, bound, post_data, post_len, 5);
    if (ZSTD_isError(compressed_size)) {
        syslog(LOG_ERR, "zts compression failed, aborting");
        free(compressed_post_data);
        free(post_data);
        return 0;
    }
    free(post_data);

    struct http_queue_entry *qe = malloc(sizeof(struct http_queue_entry));
    qe->data = compressed_post_data;
    qe->not_before_backoff = not_before;

    int push_ret = lifoq_push(sink->queue, (void*)qe, compressed_size, true, false);
    if (push_ret) {
        syslog(LOG_ERR, "HTTP Sink couldn't enqueue a %d size buffer - rejected code %d",
               post_len, push_ret);
        free(compressed_post_data);
        free(qe);
    }

    return 0;
}

static int serialize_metrics(struct http_sink* sink, metrics* m, void* data) {
    json_t** jobjects = malloc(sizeof(json_t*));
    jobjects[0] = json_object();

    const sink_config_http* httpconfig = (const sink_config_http*)sink->sink.sink_config;

    struct timeval now;
    gettimeofday(&now, NULL);

    struct cb_info info = {
        .elide = sink->elide,
        .jobjects = jobjects,
        .jobjects_count = 1,
        .config = sink->sink.global_config,
        .httpconfig = httpconfig,
        .now = now,
    };

    /* Grab the mutex state while eliding or doing other
     * operations on the shared elide map. serialize_metrics
     * can be stack invoked by statsite, especially on the
     * fast metrics path and may have more than one worker busy.
     */
    pthread_mutex_lock(&sink->sink_mutex);

    sink_elide_refresh(sink);

    /* produce a metrics json object */
    metrics_iter(m, &info, add_metrics);
    /* unlock - any shared state work should now be done until
     * lifoq
     */
    pthread_mutex_unlock(&sink->sink_mutex);

    /**
     * Compute a backoff time for this set of objects. The HTTP worker
     * will stall this long before processing queue entries, which adds a
     * local pause. Since we are using a lifoq, this means new entries
     * can establish a head of line block for requests, which is unfortunate.
     */
    struct timeval* tv = (struct timeval*) data;
    time_t not_before_backoff = 0;
    if (httpconfig->send_backoff_ms > 0) {
        double random_delay = _get_random();
        random_delay = random_delay * (double)httpconfig->send_backoff_ms;
        time_t backoff = random_delay / 1000.0;
        syslog(LOG_DEBUG, "HTTP: setting backoff time to %ld seconds", backoff);
        not_before_backoff = tv->tv_sec + (time_t)backoff;
    }

    syslog(LOG_DEBUG, "HTTP: queueing %d objects", info.jobjects_count);
    for (int i = 0; i < info.jobjects_count; i++) {
        int res = serialize_jobject(sink, info.jobjects[i], tv, not_before_backoff);
        if (res != 0) {
            syslog(LOG_NOTICE, "HTTP: couldn't package json object %d ret code %d", i, res);
            /* remember to free the remaining objects before bailing */
            for (int j = i; j < info.jobjects_count; j++) {
                json_decref(info.jobjects[j]);
            }
            free(info.jobjects);
            return res;
        }
    }
    free(info.jobjects);
    return 0;
}

/*
 * libcurl data writeback handler - buffers into a growable buffer
 */
size_t recv_cb(char *ptr, size_t size, size_t nmemb, void *userdata) {
    strbuf* buf = (strbuf*)userdata;
    /* Note: ptr is not NULL terminated, but strbuf_cat enforces a NULL */
    strbuf_cat(buf, ptr, size * nmemb);

    return size * nmemb;
}

/*
 * Attempt to check if this libcurl is using OpenSSL or NSS, which
 * differ in how ciphers are listed.
 */
static const char* curl_which_ssl(void) {
    curl_version_info_data* v = curl_version_info(CURLVERSION_NOW);
    syslog(LOG_NOTICE, "HTTP: libcurl is built with %s %s", v->version, v->ssl_version);
    if (v->ssl_version && strncmp(v->ssl_version, "NSS", 3) == 0)
        return DEFAULT_CIPHERS_NSS;
    else
        return DEFAULT_CIPHERS_OPENSSL;
}

static void http_curl_basic_setup(CURL* curl,
                                  const sink_config_http* httpconfig, struct curl_slist* headers,
                                  char* error_buffer,
                                  strbuf* recv_buf,
                                  const char* ssl_ciphers) {
    /* Setup HTTP parameters */
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, httpconfig->time_out_seconds);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, error_buffer);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, recv_buf);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, recv_cb);
    curl_easy_setopt(curl, CURLOPT_SSL_CIPHER_LIST, ssl_ciphers);
    if (headers)
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, USERAGENT);
}

/*
 * A helper to try to authenticate to an OAuth2 token endpoint
 */
static int oauth2_get_token(const sink_config_http* httpconfig, struct http_sink* sink) {
    char* error_buffer = malloc(CURL_ERROR_SIZE + 1);
    strbuf *recv_buf;
    strbuf_new(&recv_buf, 16384);

    const char* ssl_ciphers;
    if (httpconfig->ciphers)
        ssl_ciphers = httpconfig->ciphers;
    else
        ssl_ciphers = curl_which_ssl();

    CURL* curl = curl_easy_init();
    curl_easy_setopt(curl, CURLOPT_URL, httpconfig->oauth_token_url);
    http_curl_basic_setup(curl, httpconfig, NULL, error_buffer, recv_buf, ssl_ciphers);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, strlen(OAUTH2_GRANT));
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, OAUTH2_GRANT);
    curl_easy_setopt(curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
    curl_easy_setopt(curl, CURLOPT_USERNAME, httpconfig->oauth_key);
    curl_easy_setopt(curl, CURLOPT_PASSWORD, httpconfig->oauth_secret);

    CURLcode rcurl = curl_easy_perform(curl);
    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);

    int recv_len;
    char* recv_data = strbuf_get(recv_buf, &recv_len);

    if (http_code != 200 || rcurl != CURLE_OK) {
        syslog(LOG_ERR, "HTTP auth: error %d: %s %s", rcurl, error_buffer, recv_data);
        usleep(FAILURE_WAIT);
        goto exit;
    } else {
        json_error_t error;
        json_t* root = json_loadb(recv_data, recv_len, 0, &error);
        if (!root) {
            syslog(LOG_ERR, "HTTP auth: JSON load error: %s", error.text);
            goto exit;
        }
        char* token = NULL;
        if (json_unpack_ex(root, &error, 0, "{s:s}", "access_token", &token) != 0) {
            syslog(LOG_ERR, "HTTP auth: JSON unpack error: %s", error.text);
            json_decref(root);
            goto exit;
        }
        sink->oauth_bearer = strdup(token);
        json_decref(root);
        syslog(LOG_NOTICE, "HTTP auth: Got valid OAuth2 token");
    }

exit:
    curl_easy_cleanup(curl);
    free(error_buffer);
    strbuf_free(recv_buf, true);
    return 0;
}

/*
 * Let workers log who they are
 */
struct http_worker_info {
    struct http_sink* sink;
    int worker_num;
};

/*
 * A simple background worker thread which pops from the queue and tries
 * to post. If the queue is marked closed, this thread exits
 */
static void* http_worker(void* arg) {
    struct http_worker_info* info = (struct http_worker_info*)arg;
    struct http_sink* s = info->sink;
    const sink_config_http* httpconfig = (sink_config_http*)s->sink.sink_config;

    char* error_buffer = malloc(CURL_ERROR_SIZE + 1);
    strbuf *recv_buf;

    const char* ssl_ciphers;
    if (httpconfig->ciphers)
        ssl_ciphers = httpconfig->ciphers;
    else
        ssl_ciphers = curl_which_ssl();

    syslog(LOG_NOTICE, "HTTP(%d): Using cipher suite %s", info->worker_num, ssl_ciphers);

    bool should_authenticate = httpconfig->oauth_key != NULL;

    syslog(LOG_NOTICE, "HTTP(%d): Starting HTTP worker", info->worker_num);
    strbuf_new(&recv_buf, 16384);

    while(true) {
        struct http_queue_entry* queue_entry;
        void *data = NULL;
        size_t data_size = 0;
        int ret = lifoq_get(s->queue, (void**)&queue_entry, &data_size);
        data = queue_entry->data;
        if (ret == LIFOQ_CLOSED)
            goto exit;

        /* Delay sending stats until a fixed interval has elapsed */
        if (!lifoq_is_closed(s->queue)) {
            struct timeval now;
            gettimeofday(&now, NULL);
            time_t delay_for = (queue_entry->not_before_backoff - now.tv_sec);
            if (delay_for > 0)
                syslog(LOG_DEBUG, "HTTP(%d): delaying worker for %ld seconds", info->worker_num, delay_for);
            while (delay_for > 0) {
                /* Check if the queue is draining/closed, and abort sleep if needed. */
                if (lifoq_is_closed(s->queue))
                    break;
                usleep(1000000);
                delay_for -= 1;
            }
            /* Gratuitous ms level delay to jitter */
            long int ms_delay = _get_random() * 500000;
            usleep(ms_delay);
        }

        /* Hold the sink mutex for any state, such as auth cookies,
         * which may be mutated by more than one worker. */
        pthread_mutex_lock(&s->sink_mutex);
        /* Grab the current oauth bearer value so we know
           if it has changed after our request */
        const char* last_bearer = s->oauth_bearer;

        if (should_authenticate && s->oauth_bearer == NULL) {
            if (!oauth2_get_token(httpconfig, s)) {
                if (lifoq_push(s->queue, (void*)queue_entry, data_size, true, true)) {
                    syslog(LOG_ERR, "HTTP(%d): dropped data due to queue full of closed", info->worker_num);
                    if (data != NULL)
                        free(data);
                    if (queue_entry != NULL)
                        free(queue_entry);
                }
                pthread_mutex_unlock(&s->sink_mutex);
                continue;
            }
        }
        /* Build headers */
        struct curl_slist* headers = NULL;
        headers = curl_slist_append(headers, "Connection: close");

        /* Add a bearer header if needed */
        if (should_authenticate) {
            /* 30 is header preamble + fluff */
            char bearer_header[30 + strlen(s->oauth_bearer)];
            sprintf(bearer_header, "Authorization: Bearer %s", s->oauth_bearer);
            headers = curl_slist_append(headers, bearer_header);
        }

        /* Mention we are encoding this in zstd */
        headers = curl_slist_append(headers, "Content-Encoding: zstd");

        /* Release the lock after all current state has been read */
        pthread_mutex_unlock(&s->sink_mutex);

        memset(error_buffer, 0, CURL_ERROR_SIZE+1);
        CURL* curl = curl_easy_init();
        curl_easy_setopt(curl, CURLOPT_URL, httpconfig->post_url);

        http_curl_basic_setup(curl, httpconfig, headers, error_buffer, recv_buf, ssl_ciphers);

        curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, data_size);
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data);

        syslog(LOG_NOTICE, "HTTP(%d): Sending %zd bytes to %s", info->worker_num, data_size, httpconfig->post_url);
        /* Do it! */
        CURLcode rcurl = curl_easy_perform(curl);
        long http_code = 0;

        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
        if (http_code != 200 || rcurl != CURLE_OK) {

            int recv_len;
            char* recv_data = strbuf_get(recv_buf, &recv_len);

            syslog(LOG_ERR, "HTTP: error %d: %s %s", rcurl, error_buffer, recv_data);
            /* Re-enqueue data */
            if (lifoq_push(s->queue, (void*)queue_entry, data_size, true, true)) {
                if (data != NULL)
                    free(data);
                if (queue_entry != NULL)
                    free(queue_entry);
                syslog(LOG_ERR, "HTTP(%d): dropped data due to queue full of closed", info->worker_num);
            }

            /* Remove any authentication token - this will cause us to get a new one */
            pthread_mutex_lock(&s->sink_mutex);
            if (s->oauth_bearer && s->oauth_bearer == last_bearer) {
                syslog(LOG_NOTICE, "HTTP(%d): clearing Oauth bearer token", info->worker_num);
                free(s->oauth_bearer);
                s->oauth_bearer = NULL;
            }
            pthread_mutex_unlock(&s->sink_mutex);

            usleep(FAILURE_WAIT);
        } else {
            syslog(LOG_DEBUG, "HTTP(%d): success", info->worker_num);
            free(data);
            free(queue_entry);
        }

        curl_easy_cleanup(curl);
        curl_slist_free_all(headers);
        strbuf_truncate(recv_buf);
    }
exit:

    free(error_buffer);
    strbuf_free(recv_buf, true);
    return NULL;
}

static void close_sink(struct http_sink* s) {
    lifoq_close(s->queue);
    void* retval;
    for (int i = 0; i < DEFAULT_WORKERS; i++)
        pthread_join(s->worker[i], &retval);
    syslog(LOG_NOTICE, "HTTP: sink closed down with status %ld", (intptr_t)retval);
    return;
}

sink* init_http_sink(const sink_config_http* sc, const statsite_config* config) {
    struct http_sink* s = calloc(1, sizeof(struct http_sink));
    s->sink.sink_config = (const sink_config*)sc;
    s->sink.global_config = config;
    s->sink.command = (int (*)(sink*, metrics*, void*))serialize_metrics;
    s->sink.close = (void (*)(sink*))close_sink;
    s->worker = malloc(sizeof(pthread_t) * DEFAULT_WORKERS);

    unsigned int elide_generation_add = 0;
    if (rand_gather((char*)&elide_generation_add, sizeof(unsigned int)) == -1) {
        syslog(LOG_NOTICE, "HTTP: elision generation jitter not initialized");
    }

    s->elide_skip = elide_generation_add % sc->elide_interval;
    if (s->elide_skip < 0)
        s->elide_skip = 0;
    syslog(LOG_NOTICE, "HTTP: using elide skip of %d", s->elide_skip);

    sink_elide_refresh(s);

    pthread_mutex_init(&s->sink_mutex, NULL);

    syslog(LOG_NOTICE, "HTTP: using maximum queue size of %d", sc->max_buffer_size);
    lifoq_new(&s->queue, sc->max_buffer_size);
    for (int i = 0; i < DEFAULT_WORKERS; i++) {
        struct http_worker_info *info = malloc(sizeof(struct http_worker_info));
        info->sink = s;
        info->worker_num = i;
        pthread_create(&s->worker[i], NULL, http_worker, (void*)info);
    }

    return (sink*)s;
}
