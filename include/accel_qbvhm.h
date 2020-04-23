/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef QBVHM_H
#define QBVHM_H

#include <stdint.h>
#include "corona_common.h"

#ifdef __SSE__
  #include <xmmintrin.h>
#else
  #error "QBVH needs SSE! turn it on in your compile flags."
#endif

struct prims_t;

static inline void accel_print_info(FILE *fd)
{
  fprintf(fd, "accel    : qbvh with tight motion-blurred boxes and support for big scenes\n");
}

typedef union
{
  __m128 m;
  float f[4];
  unsigned int i[4];
}
qbvh_float4_t;

// double a regular qbvh node, 256 bytes
typedef struct
{
  qbvh_float4_t aabb0[6];
  qbvh_float4_t aabb1[6];
  uint64_t child[4];  // child index or, if leaf -(prim<<5|num_prims) index
  uint64_t parent;
  int64_t axis0;
  int64_t axis00;
  int64_t axis01;
}
qbvh_node_t;

typedef struct accel_t
{
  __m128 *tri_aabb;
  uint32_t num_nodes;
  uint32_t node_bufsize;
  float aabb[6];
  qbvh_node_t *tree;
  struct prims_t *prims;

  uint32_t *shadow_cache;
  uint32_t shadow_cache_last;

  uint64_t boxtests;
}
accel_t;

typedef struct
{
  int num_tris;
  int num_nodes;
  int max_depth;
  int num_leaves;
  int sum_depth;
  int num_empties;
  float parent_box;
  float SA;
}
qbvh_stats_t;

accel_t* accel_init(struct prims_t *p);
void accel_cleanup(struct accel_t *b);
void accel_build(struct accel_t *b, const char* filename);
void accel_intersect(const struct accel_t *b, const ray_t *ray, hit_t *hit);
int  accel_visible(const struct accel_t *b, const ray_t *ray);

#endif
