#include <inttypes.h>
#include <string.h>

#include "metrics.h"
#include "sink.h"
#include "utils.h"
#include "streaming.h"

struct config_time {
    const statsite_config* global_config;
    struct timeval* tv;
};

/**
 * Streaming callback to format our output
 */
static int stream_formatter(FILE *pipe, void *data, metric_type type, char *name, void *value) {
    #define STREAM(...) if (fprintf(pipe, __VA_ARGS__, (long long)tv->tv_sec) < 0) return 1;
    struct config_time* ct = data;
    struct timeval *tv = ct->tv;
    timer_hist *t;
    int i;
    char *prefix = ct->global_config->prefixes_final[type];
    switch (type) {
        case KEY_VAL:
            STREAM("%s%s|%f|%lld\n", prefix, name, *(double*)value);
            break;

        case GAUGE:
            STREAM("%s%s|%f|%lld\n", prefix, name, gauge_value(value));
            STREAM("%s%s.sum|%f|%lld\n", prefix, name, gauge_sum(value));
            STREAM("%s%s.mean|%f|%lld\n", prefix, name, gauge_mean(value));
            STREAM("%s%s.min|%f|%lld\n", prefix, name, gauge_min(value));
            STREAM("%s%s.max|%f|%lld\n", prefix, name, gauge_max(value));
            break;

        case COUNTER:
            if (ct->global_config->extended_counters) {
                STREAM("%s%s.count|%" PRIu64 "|%lld\n", prefix, name, counter_count(value));
                STREAM("%s%s.mean|%f|%lld\n", prefix, name, counter_mean(value));
                STREAM("%s%s.sum|%f|%lld\n", prefix, name, counter_sum(value));
                STREAM("%s%s.lower|%f|%lld\n", prefix, name, counter_min(value));
                STREAM("%s%s.upper|%f|%lld\n", prefix, name, counter_max(value));
                STREAM("%s%s.rate|%f|%lld\n", prefix, name, counter_sum(value) / ct->global_config->flush_interval);
            } else {
                STREAM("%s%s|%f|%lld\n", prefix, name, counter_sum(value));
            }
            break;

        case SET:
            STREAM("%s%s|%" PRIu64 "|%lld\n", prefix, name, set_size(value));
            break;

        case TIMER:
            t = (timer_hist*)value;
            STREAM("%s%s.mean|%f|%lld\n", prefix, name, timer_mean(&t->tm));
            STREAM("%s%s.lower|%f|%lld\n", prefix, name, timer_min(&t->tm));
            STREAM("%s%s.upper|%f|%lld\n", prefix, name, timer_max(&t->tm));
            STREAM("%s%s.count|%" PRIu64 "|%lld\n", prefix, name, timer_count(&t->tm));

            double quantile;
            for (i=0; i < ct->global_config->num_quantiles; i++) {
                quantile = ct->global_config->quantiles[i];
                if (quantile == 0.5) {
                    STREAM("%s%s.median|%f|%lld\n", prefix, name, timer_query(&t->tm, 0.5));
                }
                STREAM("%s%s.p%d|%f|%lld\n", prefix, name, ct->global_config->percentiles[i],
                       timer_quantile_tdigest(&t->tm, quantile));
            }
            STREAM("%s%s.rate|%f|%lld\n", prefix, name, timer_sum(&t->tm) / ct->global_config->flush_interval);
            STREAM("%s%s.sample_rate|%f|%lld\n", prefix, name, (double)timer_count(&t->tm) / ct->global_config->flush_interval);

            // Stream the histogram values
            if (t->conf) {
                STREAM("%s%s.histogram.bin_<%0.2f|%u|%lld\n", prefix, name, t->conf->min_val, t->counts[0]);
                for (i=0; i < t->conf->num_bins-2; i++) {
                    STREAM("%s%s.histogram.bin_%0.2f|%u|%lld\n", prefix, name, t->conf->min_val+(t->conf->bin_width*i), t->counts[i+1]);
                }
                STREAM("%s%s.histogram.bin_>%0.2f|%u|%lld\n", prefix, name, t->conf->max_val, t->counts[i+1]);
            }
            break;

        default:
            syslog(LOG_ERR, "Unknown metric type: %d", type);
            break;
    }
    return 0;
}

static int wrap_stream(struct sink* sink, metrics* m, void* data) {
    sink_config_stream *sc = (sink_config_stream*)sink->sink_config;
    stream_callback cb = stream_formatter;
    struct config_time ct = {
        .tv = data,
        .global_config = sink->global_config
    };
    return stream_to_command(m, &ct, cb, sc->stream_cmd);
}

/**
 * Build a stream sink - this is a basic form which has no actual parameters
 * and simply stores both configs.
 */
sink* init_stream_sink(const sink_config_stream* sc, const statsite_config* config) {
    sink* ss = calloc(1, sizeof(sink));
    ss->sink_config = (const sink_config*)sc;
    ss->global_config = config;
    ss->command = wrap_stream;
    return (sink*)ss;
}
