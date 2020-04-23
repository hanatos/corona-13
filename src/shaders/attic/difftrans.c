/*
    This file is part of corona-13.
    copyright (c) 2017 johannes hanika.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "shader.h"
#include "spectrum.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

int init(FILE *f, void **data)
{
  float *d = malloc(sizeof(float));
  *data = d;
  int i = fscanf(f, "%f", d);
  if(i < 1)
  {
    fprintf(stderr, "[difftrans] could not parse arguments! expecting: ior\n");
    d[0] = 1.5f;
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) return 1;
  return 0;
}

void cleanup(void *data)
{
  free(data);
}

float prepare(path_t *p, int v, void *data)
{ 
  float n_d = ((float *)data)[0];
  // set volume properties on next segment. we're not using this in sample(), but do it to instruct path space.
  p->v[v].interior.ior = n_d; // set fake ior
  p->v[v].material_modes = s_reflect | s_transmit | s_diffuse;
  return 1.0f;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  float wi[3];
  float wo[3];
  float n[3];
  if(e1 < e2)
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = p->e[e1].omega[k];
      wo[k] = p->e[e2].omega[k];
      n[k]  = p->v[v].hit.n[k];
    }
  }
  else
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = -p->e[e1].omega[k];
      wo[k] = -p->e[e2].omega[k];
      if(p->v[v].mode & s_transmit)
        n[k] = -p->v[v].hit.n[k];
      else
        n[k] = p->v[v].hit.n[k];
    }
  }
  const float cos_in  = -dotproduct(n, wi); // this will always be positive
  const float cos_out =  dotproduct(n, wo); // this will be positive for R, neg for T
  if(cos_in * cos_out == 0.0) return 0.0f;
  if(cos_out > 0.0f && !(p->v[v].mode & s_reflect))  return 0.0f;
  if(cos_out < 0.0f && !(p->v[v].mode & s_transmit)) return 0.0f;

  // TODO: if reflect multiply s, if transmit, multiply d
  // in projected solid angle measure
  return 0.5f / M_PI;
}


float sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.
  const float x1 = pointsampler(p, s_dim_omega_x);
  const float x2 = pointsampler(p, s_dim_omega_y);
  float s = sqrtf(x1);
  // // light tracer samples geometric normal
  // float *n = (p->v[0].mode & s_emit) ? p->v[v-1].hit.gn : p->v[v-1].hit.n;
  const float *n = p->v[v].hit.n;
  float sign = 1.0f;
  if((p->v[0].mode & s_emit) && (p->v[v].flags & s_inside))
    sign = -1.0f; // need to flip geo normal, too.
  // transmit?
  if(pointsampler(p, s_dim_scatter_mode) < .5)
  {
    p->v[v].mode = s_reflect | s_diffuse;
  }
  else
  {
    p->v[v].mode = s_transmit | s_diffuse;
    sign = - sign;
  }
  for(int k=0;k<3;k++)
    p->e[v+1].omega[k] =
      sqrtf(1.0 - x1)   * sign * n[k] +
      s*cosf(2*M_PI*x2) * p->v[v-1].hit.a[k] +
      s*sinf(2*M_PI*x2) * p->v[v-1].hit.b[k];
  p->v[v+1].pdf = 0.5/M_PI;
  return 2.0f * p->v[v].shading.rd;
}


float brdf(path_t *p, int v, void *data)
{
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);

  p->v[v].mode = s_diffuse;
  if(cos_in * cos_out < 0.0f) p->v[v].mode |= s_transmit;
  else p->v[v].mode |= s_reflect;
  return p->v[v].shading.rd / M_PI;
}

#if 0
// compute three random numbers that can be used to produce the given outgoing
// direction p->e[v+1].omega using this bsdf's sample() function.
// return 1 if this direction cannot be created at all
int inverse_sample(
    const path_t *p,
    const int v,
    float *r_omega_x,
    float *r_omega_y,
    float *r_scatter_mode,
    void *data)
{
  // TODO
}
#endif
