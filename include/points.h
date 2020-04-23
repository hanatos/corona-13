#pragma once
#include <stdint.h>
#include <stdio.h>

struct points_t;
typedef struct points_t points_t;

points_t *points_init(const uint32_t num_threads, uint64_t frame);
void points_cleanup(points_t *p);

// return uniform random in [0,1)
float points_rand(points_t *p, const uint32_t thread_num);

void points_print_info(FILE *fd);

// set the random state (for reproducible renders, if supported by the implementation)
void points_set_state(points_t *p, int thread_num, uint64_t s0, uint64_t s1);
