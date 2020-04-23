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

#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
  float mu_t;
  int mshader;
}
medium_iso_t;

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_iso_t *s = (medium_iso_t *)self->data;
  s->mshader = self - rt->shader->shader;
  return 0;
}

extern float prepare(const ray_t *ray, rayhit_t *hit, const float rr, void *data)
{ 
  medium_iso_t *s = (medium_iso_t *)data;
  hit->mu_t = s->mu_t;
  hit->mshader = s->mshader;
  return 1.0f;
}

extern float sample(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
  sample_sphere(omega_out, omega_out+1, omega_out+2, qx, qy);
  return hit->rv;
}

extern float specularity(const float *omega_in, const rayhit_t *hit, const float rr, void *data)
{
  return 0.5f;
}

extern float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return 1.0f/(4.0f*M_PI);
}

extern float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  return 1.0f/(4.0f*M_PI);
}

extern float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
  return hit->rv/(4.0f*M_PI);
}

extern int init(FILE* f, void** data)
{
  medium_iso_t *s = (medium_iso_t *)malloc(sizeof(medium_iso_t));
  *data = (void *)s;
  int i = fscanf(f, "%f", &(s->mu_t));
  if(i != 1)
  {
    fprintf(stderr, "[medium_iso] could not parse all arguments! expecting medium_iso mu_t\n");
    s->mu_t = 0.2f; // 5 dm expected path length
    return 1;
  }
  s->mu_t /= 10.0f; // adjust from 1.f->1m to 1.f->1dm
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  return 0;
}
