#pragma once

#include "pathspace.h"

typedef enum sampling_type_t
{
  st_uniform = 0,
  st_guided
} sampling_type_t;

typedef enum render_mode_t
{
  rm_learn = 0,
  rm_render
} render_mode_t;

// forward declaration of our struct:
typedef struct guided_cache_t guided_cache_t;

// initialise new cache (this will be one for all threads).
// will not cast any guide paths, just allocate storage.
guided_cache_t *guided_init(
    int num_paths);  // max number of paths in cache

// clean up memory
void guided_cleanup(guided_cache_t *c);

// clear all cdf to zero. call this before calling guided_record_path.
void guided_clear(guided_cache_t *c);

void guided_prepare_frame(
    guided_cache_t *c,
    int num_nb,
    int heap_size);

// record a path in the guided sampling cache
int guided_record_path(
    guided_cache_t *c, 
    const path_t *path, 
    const uint32_t gpath_id,
    const double measurement_contribution, 
    const double guided_pdf,
    const double uniform_pdf,
    const float uniform_ratio,
    const sampling_type_t st);

int guided_is_firefly(
    guided_cache_t *c,
    const path_t *path,
    const float throughput,
    const sampling_type_t st,
    const int normalize);

// build cdfs for sampling from what we've recorded so far.
// call this (after _clear and _record_path) before calling guided_sample).
void guided_build_cdf(
    guided_cache_t *c,
    const int num_nb,
    const float uniform_ratio,
    const render_mode_t rm);

// return pdf of whole path constructed via guided_sample() in
// projected solid angle measure (i.e. without any geometry terms).
double guided_pdf_dwp(
    const guided_cache_t *c,
    path_t *path);

// return path id to sample alongside
int guided_sample_guide_path(
    guided_cache_t *c);

uint32_t guided_num_guide_paths(
    const guided_cache_t *c);

// creates a full path from a given guide path id
// return 0 on success, 1 on failure.
int guided_sample(
    guided_cache_t *c,
    path_t *path,         // create new vertex at the end
    const uint32_t pid);  // id of guide path

void guided_debug_path_sampling(
    guided_cache_t *c,
    int num_paths);

void guided_debug_path_sampling_around_selected_pixel(
    guided_cache_t *c);

void guided_collect_stats(
    guided_cache_t* c,
    FILE* f);

void guided_gpath_hist(
    guided_cache_t* c,
    float* hist,
    int w,
    int h);

void guided_export_cache_info(
    guided_cache_t *s);
