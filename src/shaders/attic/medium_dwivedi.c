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
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

// homogeneous volume, optimised for sub surface scattering, using dwivedi sampling.
// that is, we bias the random walk to do longer steps if going back up towards the
// normal of the entry point.
// to achieve that, we present ourselves as a heterogeneous volume to the renderer.

#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>


typedef struct
{
  float g;
  float mu_t[3];
  int mshader;
}
medium_t;

const float _get_adjoint_mu_t(const path_t *p, const int e)
{
  // get adjoint entry normal:
  const float *n = 0;
  for(int v=e;v<p->length && n==0;v++)
    if(!(p->v[v].mode & s_volume))
      n = p->v[v].hit.n;
  if(!n) return p->e[e].vol.mu_t;
  const float dot = -dotproduct(n, p->e[e].omega);
  const float f = dot*0.9*p->v[e-1].interior.mu_s / p->v[e-1].interior.mu_t;
  return (1.0f-f) * p->e[e].vol.mu_t;
}

static inline float _get_mu_t(const path_t *p, const int e)
{
  // get entry normal:
  const float *n = 0;
  for(int v=e-1;v>0 && n==0;v--)
    if(!(p->v[v].mode & s_volume))
      n = p->v[v].hit.n;
  if(!n) return p->e[e].vol.mu_t;
  const float dot = dotproduct(n, p->e[e].omega);
  const float f = dot*0.9*p->v[e-1].interior.mu_s / p->v[e-1].interior.mu_t;
  return (1.0f-f) * p->e[e].vol.mu_t;
}

// switch on volume callbacks, to distinguish between hete and homo volumes:
int volume_enabled(void *data)
{
  return 1;
}

float volume_transmittance(const path_t *p, int e, void *data)
{
  return expf(-p->e[e].vol.mu_t * p->e[e].dist);
}

float volume_sample(path_t *p, int e, void *data)
{
  if(1)//p->v[e-1].interior.rv > .0f)
  { // in scattering medium
    float mu_t = _get_mu_t(p, e);
    const float rf = pointsampler(p, s_dim_free_path);
    return rt.epsilon - logf(fmaxf(rt.epsilon, 1.0f - rf))/mu_t;
  }
  return FLT_MAX;
}

float volume_pdf_adjoint(const path_t *p, int e, void *data)
{
  if(1)//p->e[e].vol.rv > .0f)
  {
    const float mu_t = _get_adjoint_mu_t(p, e);
    const float pdf = expf(-mu_t * p->e[e].dist);
    if(primid_invalid(p->v[e].hit.prim) && !(p->v[e].flags & s_environment))
    {
      return pdf * mu_t;
    }
    return pdf;
  }
  return 1.0f;
}

float volume_pdf(const path_t *p, int e, void *data)
{
  if(1)//p->e[e].vol.mu_s > .0f)
  {
    const float mu_t = _get_mu_t(p, e);
    const float pdf = expf(-mu_t * p->e[e].dist);
    if(primid_invalid(p->v[e].hit.prim) && !(p->v[e].flags & s_environment))
    {
      return pdf * mu_t;
    }
    return pdf;
  }
  return 1.0f;
}

float volume_throughput(const path_t *p, int e, void *data)
{
  if(1)//p->e[e].vol.mu_s > .0f)
  {
    // XXX this should be wrong. but dividing out the honest pdf below is a catastrophy (no idea why).
    // XXX ignoring the altered pdf here isn't all that bad, we'll be doing the same mistake the onther
    // XXX way around when going back up in the medium. of course that only holds for a semi-infinite slab..
    // return 1;
    // transmittance / pdf but bsdf in media will multiply albedo = mu_s/mu_t.
    // so we need to multiply transmittance / pdf * mu_t = expf(-mu_t d) / expf(-mu_t' d) * mu_t/mu_t'
    // and this mu_t correction factor disappears on interfaces
    const float mu_t = _get_mu_t(p, e);
    float throughput = expf(-(p->e[e].vol.mu_t - mu_t) * p->e[e].dist);
    if(primid_invalid(p->v[e].hit.prim))
      throughput *= p->e[e].vol.mu_t / mu_t;
    return throughput;
  }
  return expf(-p->e[e].vol.mu_t * p->e[e].dist);
}


int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = (medium_t *)self->data;
  s->mshader = self - rt.shader->shader;
  return 0;
}

float prepare(path_t *p, int v, void *data)
{ 
  medium_t *s = (medium_t *)data;
  // set volume properties on next segment (we pretend to be heterogeneous, so we're called on every internal vertex, too)
  const float old_mu_t = p->v[v].interior.mu_t;
  if(!(p->v[v].flags & s_environment) && primid_invalid(p->v[v].hit.prim))
  {
    p->e[v].vol.mu_t = 106.86535f * spectrum_rgb_to_p(p->lambda, s->mu_t);
    p->e[v].vol.mean_cos = s->g;
    p->v[v].material_modes = s_volume | s_glossy;
  }
  // else
  {
    p->v[v].interior.mu_t = 106.86535f * spectrum_rgb_to_p(p->lambda, s->mu_t);
    p->v[v].interior.mean_cos = s->g;
    p->v[v].interior.shader = s->mshader;
  }
  p->v[v].interior.mu_s *= p->v[v].interior.mu_t/old_mu_t;
  return 1.0f;
}

float sample(path_t *p, void *data)
{
  const int v = p->length-1;
  p->v[v].mode |= s_glossy | s_volume;
  const hit_t *hit = &p->v[v].hit;
  float out[3];
  sample_hg(p->v[v].interior.mean_cos, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), out, &p->v[v+1].pdf);
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = hit->n[k] * out[0] + hit->a[k] * out[1] + hit->b[k] * out[2];
  return p->v[v].interior.mu_s;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_volume)) return 0.0f;
  return sample_eval_hg(p->v[v].interior.mean_cos, p->e[e1].omega, p->e[e2].omega);
}

float brdf(path_t *p, int v, void *data)
{
  p->v[v].mode |= s_glossy | s_volume;
  return p->v[v].interior.mu_s * sample_eval_hg(p->v[v].interior.mean_cos, p->e[v].omega, p->e[v+1].omega);
}

int init(FILE* f, void** data)
{
  medium_t *s = (medium_t *)malloc(sizeof(medium_t));
  memset(s, 0, sizeof(medium_t));
  *data = (void *)s;
  int i = fscanf(f, "%f %f %f %f", s->mu_t, s->mu_t+1, s->mu_t+2, &s->g);
  if(i != 4)
  {
    fprintf(stderr, "[medium_dwivedi] could not parse all four arguments! expecting medium_dwivedi <mu_t rgb> <mean cosine>\n");
    for(int k=0;k<3;k++) s->mu_t[k] = 0.2f; // 2 dm expected path length
    s->g = 0.0f;
    return 1;
  }
  for(int k=0;k<3;k++)
  {
    s->mu_t[k] = 1.0f/s->mu_t[k];
    if(!isfinite(s->mu_t[k])) fprintf(stderr, "[medium_dwivedi] mean path length is zero!\n");
  }
  int res = fscanf(f, "%*[^\n]\n");

  return res == -1;
}

void cleanup(void *data)
{
  free(data);
}

