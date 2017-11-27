#ifndef TDIGEST_H
#define TDIGEST_H

#define DEFAULT_COMPRESSION 200

typedef struct centroid {
    uint64_t weight;
    double mean;
} centroid;

typedef struct point {
    double value;
    uint64_t weight;
    struct point *next;
} point;

typedef struct t_digest {
    double compression;
    int threshold;
    uint64_t size;

    uint64_t total_weight;
    double min;
    double max;

    int num_buffered_pts;
    point* buffered_pts;

    int num_centroids;
    centroid* centroids;
} t_digest;

typedef struct merging_digest {
    t_digest* t;
    centroid* centroids;
    int idx;
    double weight_so_far;
    double k1;
    double min;
    double max;
} merging_digest;

t_digest* init_tdigest(int compression);
void t_digest_add(t_digest* t, double x, uint64_t w);
double t_digest_quantile(t_digest* t, double q);
void t_digest_compress(t_digest* t);
void t_digest_free(t_digest* t);

#endif /* TDIGEST_H */