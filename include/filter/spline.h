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

#ifndef FILTER_SPLINE_H
#define FILTER_SPLINE_H

#include "corona_common.h"
#include <stdio.h>

// b = .5(t+1.5)^2         [-1.5, -.5)
//     -(t+.5)^2 + t + 1   [-.5, .5)
//     .5(t-1.5)^2          [ .5, 1.5)
// support 4x4 pixels, SIMD ?

static inline float filter_b(const float t)
{
  if(t >= -1.5 && t <= -.5)    return .5*(t+1.5f)*(t+1.5f);
  else if(t > -.5f && t <= .5) return -(t+.5)*(t+.5) + t + 1.0f;
  else if(t > .5 && t <= 1.5)  return .5*(t-1.5)*(t-1.5);
  else return 0.0f;
}

static inline void filter_accum(const float i, const float j, const float *const rad, float *pixel, const int bpp, const int offs, const int channels)
{
  int x = (int)i;
  int y = (int)j;
  float fx = i - x;
  float fy = j - y;

  const float dx = (int)(fx + .5f) - 1.5f - fx;
  const float dy = (int)(fy + .5f) - 1.5f - fy;
  float bx[4], by[4];
  for(int k=0;k<4;k++) bx[k] = filter_b(dx+k);
  for(int k=0;k<4;k++) by[k] = filter_b(dy+k);

  // printf("dx dy %f %f\n", dx, dy);
  // printf("i j %f %f x y %d %d fx fy %f %f\n", i, j, x, y, fx, fy);
  // for(int v=0;v<4;v++) for(int u=0;u<4;u++) printf("%f %f %f\n", bx[u], by[v], bx[u]*bx[v]);
  //for(float t=-2.;t<=2.;t+=0.1) printf("%f %f\n", t, filter_b(t));
  //exit(0);

  const int x0 = (int)(i + dx);
  const int y0 = (int)(j + dy);
  //printf("i j %f %f x y %d %d fx fy %f %f dx dy %f %f x0 y0 %d %d \n", i, j, x, y, fx, fy, dx, dy, x0, y0);
  const int u0 = x0 >= 0 ? 0 : - x0;
  const int v0 = y0 >= 0 ? 0 : - y0;
  const int u4 = x0 + 4 > rt.width  ? rt.width  - x0 : 4;
  const int v4 = y0 + 4 > rt.height ? rt.height - y0 : 4;
  // printf("x0 y0 %d %d u0 v0 u4 v4 %d %d %d %d\n", x0, y0, u0, v0, u4, v4);
  // if(x0 < 0 || y0 < 0 || y0 + 4 > rt.height || x0 + 4 > rt.width)
  // {
  //   float *p = pixel + 4*(x+rt.width*y) + 1;
  //   for(int k=0;k<3;k++) p[k] += rad[k];
  // }
  // else
  // {
    for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++) for(int k=0;k<channels;k++)
    {
#ifdef FILTER_ATOMIC
      common_atomic_add(pixel+bpp*(x0+u+rt.width*(y0+v)) + offs + k, rad[k]*bx[u]*by[v]);
#else
      pixel[bpp*(x0+u+rt.width*(y0+v)) + offs + k] += rad[k]*bx[u]*by[v];
#endif
    }
  //}
}

static inline void filter_print_info(FILE *fd)
{
  fprintf(fd, "filter   : spline (degree 2)");
#ifdef FILTER_ATOMIC
  fprintf(fd, ", using spinlocks");
#endif
  fprintf(fd, "\n");
}


#endif
