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

#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum.h"
#include "points.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
  float scale;
  float g;
  int mshader;
}
medium_t;

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = (medium_t *)self->data;
  s->mshader = self - rt.shader->shader;
  return 0;
}

float prepare(path_t *p, int v, void *data)
{ 
  medium_t *s = data;
  // set volume properties on next segment
  p->v[v].interior.mean_cos = s->g;
  // reverse compute mu_t from diffuse texture on the surface:
  const mf_t old_mu_t = p->v[v].interior.mu_t;
  p->v[v].interior.mu_t = mf_div(p->v[v].shading.rd, mf_set1(s->scale));
  p->v[v].interior.mu_s = mf_mul(p->v[v].interior.mu_s, mf_div(p->v[v].interior.mu_t, old_mu_t));
  if(mf_any(mf_gt(p->v[v].shading.rd, mf_set1(0.0f)))) // actually degenerate to clear medium if diffuse was 0
    p->v[v].interior.shader = s->mshader;
  // cannot do this here, we are only called upon entering the medium (homo medium),
  // would thus overwrite correct surface lobe modes!
  // p->v[v].material_modes = s_volume | s_glossy;
  return 1.0f;
}

mf_t sample(path_t *p, void *data)
{
  const int v = p->length-1;
  p->v[v].mode |= s_glossy | s_volume;
  const hit_t *hit = &p->v[v].hit;
  float out[3];
  float pdf;
  sample_hg(p->v[v].interior.mean_cos, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), out, &pdf);
  p->v[v+1].pdf = mf_set1(pdf);
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = hit->n[k] * out[0] + hit->a[k] * out[1] + hit->b[k] * out[2];
  return p->v[v].interior.mu_s;
}

int inverse_sample(
    const path_t *p,
    const int v,
    float *r_omega_x,
    float *r_omega_y,
    float *r_scatter_mode,
    void *data)
{
  float wo[3] = {
    dotproduct(p->e[v+1].omega, p->v[v].hit.n),
    dotproduct(p->e[v+1].omega, p->v[v].hit.a),
    dotproduct(p->e[v+1].omega, p->v[v].hit.b)};
  sample_inverse_hg(
      p->v[v].interior.mean_cos,
      wo,
      r_omega_x,
      r_omega_y);
  *r_scatter_mode = points_rand(rt.points, common_get_threadid());
  return 0;
}

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_volume)) return mf_set1(0.0f);
  return mf_set1(sample_eval_hg(p->v[v].interior.mean_cos, p->e[e1].omega, p->e[e2].omega));
}

mf_t brdf(path_t *p, int v, void *data)
{
  p->v[v].mode |= s_glossy | s_volume;
  return mf_mul(p->v[v].interior.mu_s, mf_set1(sample_eval_hg(p->v[v].interior.mean_cos, p->e[v].omega, p->e[v+1].omega)));
}

int init(FILE* f, void** data)
{
  medium_t *s = malloc(sizeof(medium_t));
  *data = (void *)s;
  int i = fscanf(f, "%f %f", &s->scale, &s->g);
  if(i != 2)
  {
    fprintf(stderr, "[medium poe] could not parse all arguments! expecting medium <mu_t scale> <mean cosine>\n");
    s->scale = 0.0f;
    s->g = 0.0f;
    return 1;
  }
  int res = fscanf(f, "%*[^\n]\n");
  return res == -1;
}

