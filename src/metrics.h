#ifndef METRICS_H
#define METRICS_H
#include <stdint.h>
#include "config.h"
#include "radix.h"
#include "counter.h"
#include "timer.h"
#include "gauge.h"
#include "hashmap.h"
#include "set.h"

typedef struct key_val {
    char *name;
    double val;
    struct key_val *next;
} key_val;

typedef struct {
    timer tm;

    // Support for histograms
    histogram_config *conf;
    unsigned int *counts;
} timer_hist;

typedef struct {
    hashmap *counters;  // Hashmap of name -> counter structs
    hashmap *timers;    // Map of name -> timer_hist structs
    hashmap *sets;      // Map of name -> set_t structs
    hashmap *gauges;    // Map of name -> guage struct
    key_val *kv_vals;   // Linked list of key_val structs
    double timer_eps;   // The error for timers
    int compression;    // t-digest compression
    double *quantiles;  // Array of quantiles
    uint32_t num_quants; // Size of quantiles array
    radix_tree *histograms; // Radix tree with histogram configs
    unsigned char set_precision; // The precision for sets
} metrics;

typedef int(*metric_callback)(void *data, metric_type type, char *name, void *val);

/**
 * Initializes the metrics structure.
 * @param timer_eps             the maximum error for the quantiles
 * @param timer_compression     t-digest compression
 * @param quantiles             sorted array of double quantile values, must be on (0, 1)
 * @param num_quants            number of entries in the quantiles array
 * @param histograms            radix tree with histogram settings
 * @param set_precision         the precision to use for sets
 * @param m                     metrics structure
 * @return 0 on success
 */
int init_metrics(double timer_eps, int timer_compression, double *quantiles, uint32_t num_quants, radix_tree *histograms, unsigned char set_precision, metrics *m);

/**
 * Initializes the metrics struct, with preset configurations.
 * This defaults to a epsilon of 0.01 (1% error), and quantiles at
 * 0.5, 0.95, and 0.99.
 * @return 0 on success.
 */
int init_metrics_defaults(metrics *m);

/**
 * Destroys the metrics
 * @return 0 on success.
 */
int destroy_metrics(metrics *m);

/**
 * Adds a new sampled value
 * arg type The type of the metrics
 * @arg name The name of the metric
 * @arg val The sample to add
 * @arg sample_rate The sample rate of val
 * @return 0 on success.
 */
int metrics_add_sample(metrics *m, metric_type type, char *name, double val, double sample_rate);

/**
 * Adds a value to a named set.
 * @arg name The name of the set
 * @arg value The value to add
 * @return 0 on success
 */
int metrics_set_update(metrics *m, char *name, char *value);

/**
 * Iterates through all the metrics
 * @arg m The metrics to iterate through
 * @arg data Opaque handle passed to the callback
 * @arg cb A callback function to invoke. Called with a type, name
 * and value. If the type is KEY_VAL, it is a pointer to a double,
 * for a counter, it is a pointer to a counter, and for a timer it is
 * a pointer to a timer. Return non-zero to stop iteration.
 * @return 0 on success.
 */
int metrics_iter(metrics *m, void *data, metric_callback cb);

#endif
