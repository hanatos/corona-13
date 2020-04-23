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

#pragma once
#include "corona_common.h"
#include <stdio.h>

struct prims_t;
typedef struct prims_t prims_t;
struct accel_t;
typedef struct accel_t accel_t;

// print informative description of implementation
void accel_print_info(FILE *fd);

// init new acceleration structure for given primitives (don't build yet)
accel_t* accel_init(prims_t *p);

// free memory
void accel_cleanup(accel_t *b);

// build acceleration structure, with potential file backing
void accel_build(accel_t *b, const char* filename);

// intersect ray (closest point)
void accel_intersect(const accel_t *b, const ray_t *ray, hit_t *hit);

// test visibility up to max distance
int  accel_visible(const accel_t *b, const ray_t *ray, const float max_dist);

// find closest geo intersection point to world space point at ray with parameter centre:
// ray->pos + centre * ray->dir
void accel_closest(const accel_t *b, ray_t *ray, hit_t *hit, const float centre);

// return pointer to the 6-float minxyz-maxxyz aabb
const float *accel_aabb(const accel_t *b);
