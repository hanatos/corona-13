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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
  float mu_t[81];
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
  medium_t *s = (medium_t *)data;
  // set volume properties on next segment
  p->v[v].interior.mean_cos = s->g;
  const int binf = CLAMP(81.0f * (p->lambda - 380.0f)/(780.0f-380.0f), 0, 80);
  const int bin0 = (int)binf;
  const int bin1 = CLAMP(bin0+1, 0, 80);
  const int f = binf - bin0;
  const float old_mu_t = p->v[v].interior.mu_t;
  p->v[v].interior.mu_t = s->mu_t[bin0]*(1.0f-f) + s->mu_t[bin1]*f;
  p->v[v].interior.mu_s *= p->v[v].interior.mu_t/old_mu_t;
  p->v[v].interior.shader = s->mshader;
  // cannot do this here, we are only called upon entering the medium (homo medium),
  // would thus overwrite correct surface lobe modes!
  // p->v[v].material_modes = s_volume | s_glossy;
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
  *data = (void *)s;
  char filename[512];
  int i = fscanf(f, "%s %f", filename, &s->g);
  if(i != 2)
  {
    fprintf(stderr, "[medium] could not parse all arguments! expecting medium <mu_t_file, mean cosine>\n");
    for(int k=0;k<81;k++) s->mu_t[k] = 0.2f;
    s->g = 0.0f;
    return 1;
  }
  char fname[1024];
  snprintf(fname, sizeof(fname), "%s/%s", rt.searchpath, filename);
  FILE *ff = fopen(fname, "rb");
  if(!ff)
  {
    fprintf(stderr, "[medium] could not open file `%s'!\n", filename);
    return 1;
  }
  int r = fread(s->mu_t, sizeof(float), 81, ff);
  if(r != 81)
  {
    fprintf(stderr, "[medium] could not read 81 floats from file `%s'!\n", filename);
    fclose(ff);
    return 1;
  }
  fclose(ff);
  for(int k=0;k<81;k++) s->mu_t[k] = 0.1f * s->mu_t[k];
  int res = fscanf(f, "%*[^\n]\n");
  return res == -1;
}

