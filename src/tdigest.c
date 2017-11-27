/*
 * t-digest data structure used to compute accurate percentiles, based on the MergingDigest implementation found at:
 * https://github.com/tdunning/t-digest/blob/master/src/main/java/com/tdunning/math/stats/MergingDigest.java
 */
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tdigest.h"


inline double integrated_location(double compression, double q) {
    return ((compression) * (asin(2 * (q) - 1) + M_PI / 2) / M_PI);
}

inline int float_cmp(double a, double b, double tolerance) {
    if (fabs(a - b) <= tolerance) {
        return 1;
    }
    return 0;
}

inline double max(double a, double b) {
    return a > b ? a: b;
}

inline double min(double a, double b) {
    return a < b ? a: b;
}

// utility function for qsort
static int centroid_cmp(const void *a, const void *b);
static void merge_centroid(merging_digest* args, centroid* merge);

t_digest* init_tdigest(int compression) {
    t_digest* t = (t_digest*)malloc(sizeof(t_digest));
    memset(t, 0, sizeof(t_digest));
    t->compression = compression;
    t->size = ceil(compression * M_PI / 2) + 1;
    t->threshold = 7.5 + 0.37 * compression - 2e-4 * pow(compression, 2);
    t->min = INFINITY;

    return t;
}

void t_digest_compress(t_digest* t) {
    centroid* unmerged_centroids;
    uint64_t unmerged_weight = 0;
    int i, j, num_unmerged, old_num_centroids;
    merging_digest args;

    num_unmerged = t->num_buffered_pts;
    old_num_centroids = t->num_centroids;
    if (!t->num_buffered_pts) {
        return;
    }

    unmerged_centroids = malloc(sizeof(centroid) * t->num_buffered_pts);

    for (i = 0; i < num_unmerged; i++) {
        point* p = t->buffered_pts;
        centroid* c = &unmerged_centroids[i];

        c->mean = p->value;
        c->weight = p->weight;
        unmerged_weight += c->weight;

        t->buffered_pts = p->next;
        free(p);
    }

    t->num_buffered_pts = 0;
    t->total_weight += unmerged_weight;

    qsort(unmerged_centroids, num_unmerged, sizeof(centroid), centroid_cmp);

    memset(&args, 0, sizeof(merging_digest));

    args.centroids = malloc(sizeof(centroid) * t->size);
    memset(args.centroids, 0, sizeof(centroid) * t->size);

    args.t = t;
    args.min = INFINITY;

    i = j = 0;
    while (i < num_unmerged && j < t->num_centroids) {
        centroid *a = &unmerged_centroids[i];
        centroid *b = &t->centroids[j];

        if (a->mean <= b->mean) {
            merge_centroid(&args, a);
            i++;
        } else {
            merge_centroid(&args, b);
            j++;
        }
    }

    while (i < num_unmerged) {
        merge_centroid(&args, &unmerged_centroids[i++]);
    }

    free(unmerged_centroids);

    while (j < t->num_centroids) {
        merge_centroid(&args, &t->centroids[j++]);
    }

    if (t->total_weight > 0) {
        t->min = min(t->min, args.min);

        if (args.centroids[args.idx].weight <= 0)
            args.idx--;

        t->num_centroids = args.idx + 1;
        t->max = max(t->max, args.max);
    }

    if (t->num_centroids > old_num_centroids) {
        t->centroids = realloc(t->centroids, sizeof(centroid) * t->num_centroids);
    }

    memcpy(t->centroids, args.centroids, sizeof(centroid) * t->num_centroids);

    free(args.centroids);
}

void t_digest_add(t_digest* t, double x, uint64_t w) {
    if (isnan(x) || isinf(x) || w <= 0) {
        return;
    }

    point* p = (point*)malloc(sizeof(point));
    p->value = x;
    p->weight = w;
    p->next = t->buffered_pts;

    t->buffered_pts = p;
    t->num_buffered_pts++;
    t->min = min(t->min, x);
    t->max = max(t->max, x);

    if (t->num_buffered_pts > t->threshold) {
        t_digest_compress(t);
    }
}

double t_digest_quantile(t_digest* t, double q) {
    if (t == NULL || q < 0.0 || q > 1.0) {
        return 0;
    }

    int i;
    double left, right, idx;
    centroid *a;
    centroid *b;
    centroid tmp;
    uint64_t weight_so_far;

    t_digest_compress(t);

    if (t->num_centroids == 0) {
        return NAN;
    }

    if (t->num_centroids == 1)
        return t->centroids[0].mean;

    if (float_cmp(q, 0.0, FLT_EPSILON)) {
        return t->min;
    }

    if (float_cmp(q, 1.0, FLT_EPSILON)) {
        return t->max;
    }

    idx = q * t->total_weight;

    weight_so_far = 0;
    b = &tmp;
    b->mean = t->min;
    b->weight = 0;
    right = t->min;

    for (i = 0; i < t->num_centroids; i++) {
        centroid* c = &t->centroids[i];
        a = b;
        left = right;

        b = c;
        right = (b->weight * a->mean + a->weight * b->mean) / (a->weight + b->weight);

        if (idx < weight_so_far + a->weight) {
            double p = (idx - weight_so_far) / a->weight;
            return left * (1 - p) + right * p;
        }

        weight_so_far += a->weight;
    }

    left = right;
    a = b;
    right = t->max;

    if (idx < weight_so_far + a->weight) {
        double p = (idx - weight_so_far) / a->weight;
        return left * (1 - p) + right * p;
    }

    return t->max;
}

void t_digest_free(t_digest* t) {
    while (t->buffered_pts) {
        point* p = t->buffered_pts;
        t->buffered_pts = t->buffered_pts->next;
        free(p);
    }

    if (t->centroids) {
        free(t->centroids);
    }

    free(t);
}

static int centroid_cmp(const void *a, const void *b) {
    centroid* c1 = (centroid*) a;
    centroid* c2 = (centroid*) b;

    return (c1->mean == c2->mean) ? 0 : (c1->mean < c2->mean) ? -1 : 1;
}

static void merge_centroid(merging_digest* args, centroid* merge) {
    double k2;
    centroid* c = &args->centroids[args->idx];

    args->weight_so_far += merge->weight;
    k2 = integrated_location(args->t->compression, args->weight_so_far / args->t->total_weight);

    if (k2 - args->k1 > 1 && c->weight > 0) {
        args->idx++;
        args->k1 = integrated_location(args->t->compression,
                                       (args->weight_so_far - merge->weight) / args->t->total_weight);
    }

    c = &args->centroids[args->idx];
    c->weight += merge->weight;
    c->mean += (merge->mean - c->mean) * merge->weight / c->weight;

    if (merge->weight > 0) {
        args->min = min(merge->mean, args->min);
        args->max = max(merge->mean, args->max);
    }
}
