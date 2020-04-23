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

#ifndef FILTER_BILIN_H
#define FILTER_BILIN_H

#include "corona_common.h"

static inline void filter_accum(const float i, const float j, const float *const rad, float *pixel, const int bpp, const int offs, const int channels)
{
  if((i < rt.width-1) && (j < rt.height-1))
  {
    int x = (int)i;
    int y = (int)j;
    float fx = i - x;
    float fy = j - y;
    float *pixel00 = pixel + bpp*(x+rt.width*y)         + offs;
    float *pixel01 = pixel + bpp*((x+1)+rt.width*y)     + offs;
    float *pixel10 = pixel + bpp*(x+rt.width*(y+1))     + offs;
    float *pixel11 = pixel + bpp*((x+1)+rt.width*(y+1)) + offs;
    for(int k=0;k<channels;k++)
    {
#ifdef FILTER_ATOMIC
      common_atomic_add(pixel00+k, (1-fx)*(1-fy)*rad[k]);
      common_atomic_add(pixel01+k, (fx)*(1-fy) * rad[k]);
      common_atomic_add(pixel10+k, (1-fx)*(fy) * rad[k]);
      common_atomic_add(pixel11+k, (fx)*(fy)   * rad[k]);
#else
      pixel00[k] += (1-fx)*(1-fy)*rad[k];
      pixel01[k] += (fx)*(1-fy) * rad[k];
      pixel10[k] += (1-fx)*(fy) * rad[k];
      pixel11[k] += (fx)*(fy)   * rad[k];
#endif
    }
  }
  else
  {
    float *p = pixel + bpp*((int)i+rt.width*(int)j) + offs;
    for(int k=0;k<channels;k++)
    {
#ifdef FILTER_ATOMIC
      common_atomic_add(p+k, rad[k]);
#else
      p[k] += rad[k];
#endif
    }
  }
}

static inline void filter_print_info(FILE *fd)
{
  fprintf(fd, "filter   : bilinear interpolation");
#ifdef FILTER_ATOMIC
  fprintf(fd, ", using spinlocks");
#endif
  fprintf(fd, "\n");
}

#endif
