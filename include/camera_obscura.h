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
#ifndef CAMERA_H
#define CAMERA_H

#include "camera_common.h"

#include <math.h>

static inline void camera_print_info(FILE *f)
{
  fprintf(f, "camera   : obscura (pin-hole)\n");
  fprintf(f, "  pos    : %f %f %f\n", rt.cam->pos[0], rt.cam->pos[1], rt.cam->pos[2]);
  float dir[3] = {0, 0, 1.0};
  quaternion_transform(&(rt.cam->orient), dir);
  fprintf(f, "  dir    : %f %f %f\n", dir[0], dir[1], dir[2]);
  fprintf(f, "  quat   : %f, (%f %f %f)\n", rt.cam->orient.w, rt.cam->orient.x[0], rt.cam->orient.x[1], rt.cam->orient.x[2]);
  fprintf(f, "  crop   : %.1f\n", rt.cam->crop_factor);
  fprintf(f, "  film   : %dmm x %dmm\n", (int)(rt.cam->film_width*100.0), (int)(rt.cam->film_height*100.0 + 0.5f));
  fprintf(f, "  f-len  : %dmm\n", (int)(rt.cam->focal_length*100.0 + 0.5f));
}

//static inline void camera_measure(const camera_t *c, const float x, const float y, float *E)
static inline void camera_measure(const camera_t *c, const ray_t *ray, float *E)
{
  const float dot = dotproduct(ray->dir, c->dir);
  const float f = dot*dot*dot*dot/(c->focal_length*c->focal_length);
  for(int k=0;k<3;k++) E[k] = f*c->film_width*c->film_height;
}

static inline int camera_get_pixel(const camera_t *c, const rayhit_t *hit, float *i, float *j, rayhit_t *lensp)
{
  // project onto basis vectors.
  float v[3];
  for(int k=0;k<3;k++) v[k] = hit->hit[k] - c->pos[k];
  const float s = c->focus/dotproduct(c->dir, v);
  for(int k=0;k<3;k++) v[k] = v[k]*s - c->f_ul[k];
  //printf("unnormalised dir        : %f %f %f\n", v[0], v[1], v[2]);
  const float inv_fr_len = (c->width*c->focal_length)/(c->film_width * c->focus);
  const float inv_fd_len = (c->height*c->focal_length)/(c->film_height * c->focus);
  *i = dotproduct(v, c->f_r)*inv_fr_len*inv_fr_len;
  *j = dotproduct(v, c->f_d)*inv_fd_len*inv_fd_len;
  // bounds - 1 for filters using (int)i and (int)i + 1
  if(*i < 0 || *j < 0 || *i >= c->width || *j >= c->height) return 1;
  // construct lens hit:
  for(int k=0;k<3;k++)
  {
    lensp->normal[k] = c->dir[k];
    lensp->hit[k] = c->pos[k];
  }
  return 0;
}

static inline void camera_get_ray(camera_t *c, const float i, const float j, ray_t *ray, rayhit_t *hit)
{
  //create vector:
  for(int k=0;k<3;k++) ray->dir[k] = c->f_ul[k] + i*c->f_r[k] + j*c->f_d[k];

  float len = sqrtf(ray->dir[0]*ray->dir[0] + ray->dir[1]*ray->dir[1] + ray->dir[2]*ray->dir[2]);
  for(int k=0;k<3;k++) ray->dir[k] /= len;

  for(int k=0;k<3;k++) ray->pos[k] = c->pos[k];
  hit->prim = -1;
  hit->inside = 0;
  hit->ior = 1.0f;
  hit->mu_t = 0.0f;
  for(int k=0;k<3;k++) ray->invdir[k] = 1.0f/ray->dir[k];
}

#endif
