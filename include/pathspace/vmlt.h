#pragma once

#include "pathspace.h"

// access to brain of vmlt for mutation strategies

typedef struct vmlt_thr_t
{
  uint64_t num_rejects;
  float acceptance;
  float accum;
  int large_step;
  int mutation;
  uint64_t num_mutations;
  uint64_t accepted[20];
  uint64_t rejected[20];
  path_t *curr_path; // remember last current path so we can splat it in finalize()
}
vmlt_thr_t;

// modular struct to store function pointers to mutation strategies
typedef void* (*vmlt_init_t)();
typedef void (*vmlt_cleanup_t)(void *);
typedef float (*vmlt_suitability_t)(const path_t *p, void *);
typedef float (*vmlt_mutate_t)(path_t *curr, path_t *tent, void *);
typedef void (*vmlt_print_info_t)(FILE *, void *);
typedef struct vmlt_mutation_t
{
  vmlt_init_t init;
  vmlt_cleanup_t cleanup;
  vmlt_suitability_t suitability;
  vmlt_mutate_t mutate;
  vmlt_print_info_t print_info;
  void *data;
}
vmlt_mutation_t;

typedef struct vmlt_t
{
  vmlt_thr_t *t;
  int mutations;
  uint64_t max_num_rejects;
  vmlt_mutation_t mutation[20]; // max mutation strategies is 20
}
vmlt_t;

void vmlt_register(
    vmlt_t *p,
    vmlt_init_t init,
    vmlt_cleanup_t cleanup,
    vmlt_suitability_t s,
    vmlt_mutate_t m,
    vmlt_print_info_t i);

void vmlt_mutate(vmlt_t *s, path_t *curr, path_t *tent);
int vmlt_accept(vmlt_t *s, path_t *curr, path_t *tent);
vmlt_t *vmlt_init();
void vmlt_cleanup(vmlt_t *s);
