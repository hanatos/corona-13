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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "shader.h"
#include "accel.h"
#include "rgb2spec.h"
#include "spectrum.h"
#include "sampler_common.h"

typedef struct
{
  float coeff[3];
  float scale;
}
const_t;

void cleanup(void *data)
{
  free(data);
}

mf_t eval(path_t *p, int v, void *data)
{
  const_t *t = data;
  return mf_mul(mf_rgb2spec(t->coeff, p->lambda), mf_set1(t->scale));
}

mf_t sample(path_t *p, void *data)
{
  const_t *t = data;
  float x1, x2;
  if(p->v[0].mode & s_emit)
  { // light ray
    x1 = pointsampler(p, s_dim_edf_x);
    x2 = pointsampler(p, s_dim_edf_y);
  }
  else
  { // next event
    x1 = pointsampler(p, s_dim_nee_x);
    x2 = pointsampler(p, s_dim_nee_y);
  }
  const int v = p->length;
  sample_sphere(p->e[v].omega, p->e[v].omega+1, p->e[v].omega+2, x1, x2);
  p->e[v].dist = FLT_MAX;
  p->v[v].shading.roughness = 1.0f;
  p->v[v].flags = s_environment;
  p->v[v].mode = s_emit;
  p->v[v].pdf = mf_set1(1.0f/(4.0f*M_PI));
  p->v[v].shading.em = mf_mul(mf_rgb2spec(t->coeff, p->lambda), mf_set1(t->scale));
  if(p->length)
  {
    const float *aabb = accel_aabb(rt.accel);
    const float far = aabb[3] + aabb[4] + aabb[5]
      - aabb[0] - aabb[1] - aabb[2];
    for(int k=0;k<3;k++)
      p->v[v].hit.x[k] = p->v[v-1].hit.x[k] + far*p->e[v].omega[k];
  }
  for(int k=0;k<3;k++)
    p->v[v].hit.n[k] = p->v[v].hit.gn[k] = - p->e[v].omega[k];
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = -1;
  return mf_div(p->v[v].shading.em, p->v[v].pdf);
}

mf_t pdf(path_t *p, int v, void *data)
{
  return mf_set1(1.0f/(4.0f*M_PI));
}

int init(FILE* s, void **data)
{
  const_t *t = malloc(sizeof(*t));
  *data = t;
  char line[1024];
  int dreggn = 0;
  dreggn = fscanf(s, "%[^\n]\n", line);
  float col[3] = {1.0f, 1.0f, 1.0f};
  t->scale = 1.0f;
  dreggn += sscanf(line, "%f %f %f %f", col, col+1, col+2, &t->scale);
  t->scale *= spectrum_rgb_to_coeff(col, t->coeff);
  return dreggn == -1;
}

