#pragma once

#include "mf.h"
#include <stdint.h>
#include <stdio.h>

struct render_t;
typedef struct render_t render_t;
struct render_tls_t;
typedef struct render_tls_t render_tls_t;
struct path_t;

struct render_t *render_init();
void render_cleanup(render_t *r);
struct render_tls_t *render_tls_init();
void render_tls_cleanup(render_tls_t *r);

// print descriptive string
void render_print_info(FILE *fd);

// sample path with given index
void render_sample_path(uint64_t index);

// clear frame
void render_clear();

// splat a path with given contribution
void render_splat(const struct path_t *p, const mf_t value);
