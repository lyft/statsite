#include <check.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "tdigest.h"

#define EPSILON 21474836

static int check_answer(double computed, double answer, double tol) {
    if (isnan(computed) || isnan(answer))
        return 1;

    if (computed < answer + tol && computed > answer - tol) {
        return 1;
    } else {
        fprintf(stderr, "[ERROR] expected %.9g, got %.9g\n", answer, computed);
        return 0;
    }
}

START_TEST(test_td_init_and_destroy)
    {
        t_digest* td;
        td = init_tdigest(DEFAULT_COMPRESSION);
        fail_unless(td != NULL);
        t_digest_free(td);
    }
END_TEST

START_TEST(test_td_init_add_destroy)
    {
        t_digest* td;
        td = init_tdigest(DEFAULT_COMPRESSION);
        fail_unless(td != NULL);

        t_digest_add(td, 100.0, 1.0);
        fail_unless(t_digest_quantile(td, 1.0) == 100);
        fail_unless(t_digest_quantile(td, 0.0) == 100);
        t_digest_free(td);
    }
END_TEST

START_TEST(test_td_stream_add_check_upper_lower)
    {
        int i, num = 7;
        double data_stream[7] = {5.4, -3.2, 2.0, 2.5, 2.7, 9.8, -5.6};
        double data_stream_sorted[7] = {-5.6, -3.2, 2.0, 2.5, 2.7, 5.4, 9.8};

        t_digest* td = init_tdigest(DEFAULT_COMPRESSION);
        for (i = 0; i < num; i++) {
            t_digest_add(td, data_stream[i], 1);
        }
        fail_unless(check_answer(t_digest_quantile(td, 1.0), data_stream_sorted[6], EPSILON), "Incorrect upper");
        fail_unless(check_answer(t_digest_quantile(td, 0.0), data_stream_sorted[0], EPSILON), "Incorrect lower");
    }
END_TEST

START_TEST(test_td_init_add_loop_random_query_destroy)
    {
        t_digest* td;
        td = init_tdigest(DEFAULT_COMPRESSION);
        fail_if(td == NULL);

        srandom(42);
        for (int i = 0; i < 100000; i++) {
            t_digest_add(td, random(), 1.0);
        }

        double val = t_digest_quantile(td, 0.5);
        fail_unless(check_answer(val, 1073741823, EPSILON), "Incorrect mean");

        val = t_digest_quantile(td, 0.90);
        fail_unless(check_answer(val, 1932735282, EPSILON), "Incorrect p90");

        val = t_digest_quantile(td, 0.99);
        fail_unless(check_answer(val, 2126008810, EPSILON), "Incorrect p99");

        t_digest_free(td);
    }
END_TEST


START_TEST(test_td_init_add_loop_rev_query_destroy)
    {
        t_digest* td;
        td = init_tdigest(DEFAULT_COMPRESSION);
        fail_if(td == NULL);

        for (int i = 100000; i > 0; i--) {
            t_digest_add(td, i, 1.0);
        }

        double val = t_digest_quantile(td, 0.5);
        fail_unless(check_answer(val, 50000.0, 1000), "Incorrect mean");

        val = t_digest_quantile(td, 0.9);
        fail_unless(check_answer(val, 90000.0, 1000), "Incorrect p90");

        val = t_digest_quantile(td, 0.99);
        fail_unless(check_answer(val, 99000.0, 1000), "Incorrect p99");

        t_digest_free(td);
    }
END_TEST

START_TEST(test_td_init_add_negative_query_destroy)
    {
        t_digest* td;
        td = init_tdigest(DEFAULT_COMPRESSION);
        fail_if(td == NULL);

        t_digest_add(td, -100.0, 1.0);

        double val = t_digest_quantile(td, 0.5);
        fail_unless(val == -100.0);

        t_digest_free(td);
    }
END_TEST
