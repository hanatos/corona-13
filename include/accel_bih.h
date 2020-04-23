/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef BIH_H
#define BIH_H

#include "prims.h"


static inline void accel_print_info(FILE *fd)
{
  fprintf(fd, "accel    : bih\n");
}

typedef struct
{
  float clip[2];
  unsigned int data;
}
bih_node_t;

typedef struct accel_t
{
  unsigned int num_nodes;
  unsigned int num_clipnodes;
  unsigned int node_bufsize;
  float aabb[6];
  struct prims_t *prims;
  bih_node_t *tree;
  long long boxtests;
}
accel_t;

#if 0
static inline char accel_is_leaf(bih_node_t* n)
{
  return (n->data & 3) == 3;
}

static inline char accel_is_clipnode(bih_node_t* n)
{
  return (n->data & 4);
}

static inline bih_node_t* accel_children(bih_node_t* n, accel_t *t)
{
  return t->tree + (n->data >> 3);
}

static inline unsigned int accel_dim(bih_node_t* n)
{
  return n->data & 3;
}

static inline unsigned int accel_prims(bih_node_t *n)
{
  return n->data >> 3;
}

static inline unsigned int accel_num(bih_node_t *n)
{
  return *(unsigned int*)n;
}
#endif

accel_t* accel_init(prims_t *p);
void accel_cleanup(struct accel_t *b);
void accel_build(struct accel_t *b, const char* filename);
void accel_intersect(const struct accel_t *b, const ray_t *ray, hit_t *hit);
int  accel_visible(const struct accel_t *b, const ray_t *ray);

#endif
