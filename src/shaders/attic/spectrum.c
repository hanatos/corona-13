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
#include "sampler_common.h"
#include "shader.h"
#include "spectrum.h"

const int SLOT_DIFFUSE = 0;
const int SLOT_SPECULAR = 1;
const int SLOT_EMISSION = 2;
const int SLOT_VOLUME = 3;

typedef struct
{
	int lambda_min, lambda_step, lambda_num;
  float *spectrum, *cdf;
  float roughness, max;
  int slot;
}
spectrum_t;

static inline float get_spectrum(const spectrum_t *s, float lambda)
{
  int l = (lambda - s->lambda_min)/(float)s->lambda_step;
  if(l < 0 || l >= s->lambda_num) return 0.0f;
  return s->spectrum[l];
}

int init(FILE *s, void **data)
{
  spectrum_t *t = (spectrum_t *)malloc(sizeof(spectrum_t));
  *data = t;
  //// read sample shader number
  char c;
  char filename[512];
  t->slot = SLOT_DIFFUSE;
  t->roughness = 1.0f; // cos^1
  float scale = 1.0f;
  if(fscanf(s, " %c %s %f %f", &c, filename, &scale, &(t->roughness)) < 3)
  {
    fprintf(stderr, "[spectrum] shader could not read all parameters! expecting: <d|s|e|v> <spectrum.txt> <scale> [roughness]\n");
    return 1;
  }
  if(c == 's') t->slot = SLOT_SPECULAR;
  else if(c == 'e') t->slot = SLOT_EMISSION;
  else if(c == 'v') t->slot = SLOT_VOLUME;
  FILE *fd = fopen(filename, "rd");
  if(!fd)
  {
    printf("[spectrum] ERROR: could not open %s!\n", filename);
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(fd, "%d %d %d\n", &(t->lambda_min), &(t->lambda_step), &(t->lambda_num));
  t->spectrum = (float *)malloc(sizeof(float)*t->lambda_num);
  t->cdf = (float *)malloc(sizeof(float)*t->lambda_num);
  // printf("[spectrum] scale %f [%d %d %d]\n", scale, t->lambda_min, t->lambda_step, t->lambda_num);
  int k = 0;
  t->max = 0.0f;
  while (!feof(fd) && k < t->lambda_num)
  {
    dreggn = fscanf(fd, "%f", t->spectrum + k++);
    t->spectrum[k-1] *= scale;
    t->max = fmaxf(t->spectrum[k-1], t->max);
    dreggn = fscanf(fd, "%*[^\n]\n");
  }
  fclose(fd);

  dreggn = fscanf(s, "%*[^\n]\n");
  return dreggn == -1;
}

void cleanup(void *data)
{
  spectrum_t *s = (spectrum_t *)data;
  free(s->spectrum);
  free(s->cdf);
  free(s);
}

float emission(float *roughness, float *em, const float lambda, void *data)
{
  spectrum_t *s = (spectrum_t *)data;
  if(s->slot == SLOT_EMISSION)
  {
    *roughness = s->roughness;
    *em = get_spectrum(s, lambda);
    return s->max;
  }
  else return 0.0f;
}

float sample_emission(float *lambda, float *roughness, const float rand, void *data)
{
  spectrum_t *s = (spectrum_t *)data;
  int l = sample_cdf(s->cdf, s->lambda_num, rand);
  float em, jitter = rand*1337.0;
  jitter -= (int)jitter;
  *lambda = (l+jitter)*s->lambda_step + s->lambda_min;
  emission(roughness, &em, *lambda, data);
  return em;
}

float prepare(path_t *p, int v, void *data)
{
  spectrum_t *s = (spectrum_t *)data;
  p->v[v].shading.roughness = s->roughness;
  float power = get_spectrum(s, p->lambda);
  if(s->slot == SLOT_DIFFUSE) p->v[v].shading.rd = power;
  else if(s->slot == SLOT_SPECULAR) p->v[v].shading.rs = power;
  else if(s->slot == SLOT_VOLUME)
  {
    p->v[v].interior.mu_s = power;
    p->v[v].interior.mu_t = 1.0f;
  }
  else p->v[v].shading.em = power;//if(s->slot == SLOT_EMISSION) hit->em = power;
  return 1.0f;
}

