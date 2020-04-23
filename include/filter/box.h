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

#pragma once

#include "corona_common.h"
#include "framebuffer.h"

static inline void filter_box_splat(
    framebuffer_t *fb,
    const float i,
    const float j,
    const float *const col)
{
  float *p = fb->fb + fb->header->channels*((int)i+fb->header->width*(int)j);
  for(int k=0;k<fb->header->channels;k++) 
#ifdef FILTER_ATOMIC
    common_atomic_add(p+k, col[k]);
#else
    p[k] += col[k];
#endif
}

static inline void filter_box_splat4(
    framebuffer_t *fb,
    const float i,
    const float j,
    const __m128 col4)
{
  common_atomic_add128(
      fb->fb + fb->header->channels*((int)i+fb->header->width*(int)j),
      col4);
}

static inline void filter_box_print_info(FILE *fd)
{
  fprintf(fd, "filter   : box");
#ifdef FILTER_ATOMIC
  fprintf(fd, ", using spinlocks");
#endif
  fprintf(fd, "\n");
}
