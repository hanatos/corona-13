/*
    This file is part of corona-13.

    copyright (c) 2018 johannes hanika.

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

#ifndef FILTER_BLACKMAN_HARRIS_H
#define FILTER_BLACKMAN_HARRIS_H

#include "corona_common.h"
#include <stdio.h>

// support 4x4 pixels

static inline float filter_bh_w(const float n)
{
  const float NN = 4.0f;
  if(n > NN-1.0f || n < 0.0f) return 0.0f;
  const float a0 = 0.35875;
  const float a1 = 0.48829;
  const float a2 = 0.14128;
  const float a3 = 0.01168;
  const float N_1 = 1.0f/(NN-1.0f);
  const float cos1 = cosf(2.0f*M_PI*n*N_1);
  const float cos2 = cosf(4.0f*M_PI*n*N_1);
  const float cos3 = cosf(6.0f*M_PI*n*N_1);
  return a0 - a1*cos1 + a2*cos2 - a3*cos3;
}

static inline void filter_blackmanharris_splat(
    framebuffer_t *fb,
    const float i,
    const float j,
    const float *const col)
{
  const int wd = fb->header->width, ht = fb->header->height;
  const int x0 = (int)(i - 1.5f);
  const int y0 = (int)(j - 1.5f);
  const int u0 = x0 >= 0 ? 0 : - x0;
  const int v0 = y0 >= 0 ? 0 : - y0;
  const int u4 = x0 + 4 > wd ? wd - x0 : 4;
  const int v4 = y0 + 4 > ht ? ht - y0 : 4;

  float weight = 0.0f;
  for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++)
  {
    const float uu = (x0 + u + .5f) - i, vv = (y0 + v + .5f) - j;
    const float r = sqrtf(uu*uu + vv*vv);
    const float f = filter_bh_w(r+1.5f);
    weight += f;
  }
  if(weight <= 0) return;
  weight = 1.0f/weight;
  for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++)
  {
    const float uu = (x0 + u + .5f) - i, vv = (y0 + v + .5f) - j;
    const float r = sqrtf(uu*uu + vv*vv);
    const float f = weight * filter_bh_w(r+1.5f);
    const float col2[3] = {col[0]*f, col[1]*f, col[2]*f};
    filter_box_splat(fb, x0+u, y0+v, col2);
  }
}


static inline void filter_blackmanharris_splat4(
    framebuffer_t *fb,
    const float i,
    const float j,
    const __m128 col4)
{
  const int wd = fb->header->width, ht = fb->header->height;
  const int x0 = (int)(i - 1.5f);
  const int y0 = (int)(j - 1.5f);
  const int u0 = x0 >= 0 ? 0 : - x0;
  const int v0 = y0 >= 0 ? 0 : - y0;
  const int u4 = x0 + 4 > wd ? wd - x0 : 4;
  const int v4 = y0 + 4 > ht ? ht - y0 : 4;

  float weight = 0.0f;
  for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++)
  {
    const float uu = (x0 + u + .5f) - i, vv = (y0 + v + .5f) - j;
    const float r = sqrtf(uu*uu + vv*vv);
    const float f = filter_bh_w(r+1.5f);
    weight += f;
  }
  if(weight <= 0) return;
  weight = 1.0f/weight;
  __m128 rad4 = _mm_mul_ps(_mm_set1_ps(weight), col4);
  for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++)
  {
    const float uu = (x0 + u + .5f) - i, vv = (y0 + v + .5f) - j;
    const float r = sqrtf(uu*uu + vv*vv);
    const float f = filter_bh_w(r+1.5f);
    __m128 f4 = _mm_set1_ps(f);
    __m128 radf = _mm_mul_ps(rad4, f4);
    common_atomic_add128(fb->fb+4*(x0+u+wd*(y0+v)), radf);
  }
}

static inline void filter_blackmanharris_print_info(FILE *fd)
{
  fprintf(fd, "filter   : blackman harris");
#ifdef FILTER_ATOMIC
  fprintf(fd, ", using spinlocks");
#endif
  fprintf(fd, "\n");
}

#endif
