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
#include "spectrum.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct
{
  // absorbtion in dye
  float *mu_t;
  // photo luminiscence spectrum, cdf
  float *pl_cdf;
  int pl_min;
  int pl_step;
  int pl_num;
  int mu_min;
  int mu_step;
  int mu_num;
  int mshader;
  const rt_t *rt;
}
medium_fluko_t;

float get_spectrum(float *s, int min, int step, int num, float lambda)
{
  int l = (lambda - min)/step;
  if(l < 0 || l >= num) return 0.0f;
  return s[l];
}

float get_spectrum2(float *s, int min, int step, int num, float lambda)
{
  int l = (lambda - min)/step;
  if(l < 0 || l >= num) return 0.0f;
  if(l == 0) return s[l];
  return (s[l] - s[l-1])/step;
}

#define QUANTUM_EFFICIENCY 0.75f// 0.95f

int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_fluko_t *s = (medium_fluko_t *)self->data;
  s->rt = rt;
  s->mshader = self - rt->shader->shader;
  return 0;
}

extern float prepare(const ray_t *ray, rayhit_t *hit, const float rr, void *data)
{ 
  medium_fluko_t *s = (medium_fluko_t *)data;
  hit->mu_t = get_spectrum(s->mu_t, s->mu_min, s->mu_step, s->mu_num, hit->lambda);
  hit->mshader = s->mshader;
  return 1.0f;
}

extern float sample_adj(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
  medium_fluko_t *s = (medium_fluko_t *)data;
  sample_sphere(omega_out, omega_out+1, omega_out+2, qx, qy);

  // stokes shift:
  // adjoint (pt)
  const float pl = get_spectrum2(s->pl_cdf, s->pl_min, s->pl_step, s->pl_num, hit->lambda);
  float pdf;
  hit->lambda = spectrum_sample_lambda(rr, &pdf);
  // assert(hit->lambda >= 400 && hit->lambda <= 700);
  // printf("sample adj pl: %f\n", pl*QUANTUM_EFFICIENCY*(700.0f-400.0f));
  return pl*QUANTUM_EFFICIENCY/pdf;
}

extern float sample(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
  medium_fluko_t *s = (medium_fluko_t *)data;
  sample_sphere(omega_out, omega_out+1, omega_out+2, qx, qy);

  float pdf;
  hit->lambda = spectrum_sample_lambda(rr, &pdf);
  const float pl = get_spectrum2(s->pl_cdf, s->pl_min, s->pl_step, s->pl_num, hit->lambda);

  // printf("sample     pl: %f\n", pl);
  // assert(hit->lambda >= 400 && hit->lambda <= 700);
  return pl*QUANTUM_EFFICIENCY/pdf;
  // stokes shift:
  int k = sample_cdf(s->pl_cdf, s->pl_num, rr);
  hit->lambda = s->pl_min + s->pl_step*k;
  // if(hit->lambda < 400 || hit->lambda > 700) printf("lmabd %f\n", hit->lambda);
  // assert(hit->lambda >= 400 && hit->lambda <= 700);
  return QUANTUM_EFFICIENCY;//*1.0f/(s->pl_cdf[k] - (k > 0 ? s->pl_cdf[k-1] : 0.0f));
}

extern float specularity(const float *omega_in, const rayhit_t *hit, const float rr, void *data)
{
  return 0.5f;
}

extern float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return 0.0f;
}

extern float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  return 0.0f;
  // return 1.0f/(4.0f*M_PI); // * lambda pdf?
}

extern float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
  return 0.0f; // prob to sample same lambda is zero :(
}

FILE *fopen_rt(const char *filename, const char *mode, const rt_t *rt)
{
  FILE* f = fopen(filename, mode);
  if(!f)
  {
    char name[512];
    sprintf(name, "%s/%s", rt->searchpath, filename);
    f = fopen(name, "rb");
  }
  return f;
}

extern void cleanup(void *data)
{
  medium_fluko_t *s = (medium_fluko_t *)data;
  free(s->mu_t);
  free(s->pl_cdf);
  free(s);
}

extern int init(FILE* fd, void** data)
{
  medium_fluko_t *s = (medium_fluko_t *)malloc(sizeof(medium_fluko_t));
  *data = (void *)s;
  char sigmaf[512], plf[512];
  int i = fscanf(fd, "%s %s", sigmaf, plf);
  if(i != 2)
  {
    fprintf(stderr, "[medium_fluko] could not parse all arguments! expecting mu_t.txt pl.txt\n");
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(fd, "%*[^\n]\n");
  FILE *f = fopen(sigmaf, "rb");//, s->rt);
  if(!f) 
  {
    fprintf(stderr, "[medium_fluko] could not read %s\n", sigmaf);
    return 1;
  }
  dreggn = fscanf(f, "%d %d %d\n", &(s->mu_min), &(s->mu_step), &(s->mu_num));
  s->mu_num -= s->mu_min;
  s->mu_num /= s->mu_step;
  s->mu_num += 1;
  s->mu_t = (float *)malloc(sizeof(float)*s->mu_num);
  for(int k=0;k<s->mu_num;k++)
  {
    dreggn = fscanf(f, "%f\n", s->mu_t + k);
    s->mu_t[k] *= 100.0f;
  }
  fclose(f);
  f = fopen(plf, "rb");//, s->rt);
  if(!f) 
  {
    fprintf(stderr, "[medium_fluko] could not read %s\n", plf);
    return 1;
  }
  dreggn = fscanf(f, "%d %d %d\n", &(s->pl_min), &(s->pl_step), &(s->pl_num));
  s->pl_num -= s->pl_min;
  s->pl_num /= s->pl_step;
  s->pl_num += 1;
  s->pl_cdf = (float *)malloc(sizeof(float)*s->pl_num);
  // const int smooth = 20;
  for(int k=0;k<s->pl_num;k++)
  {
    dreggn = fscanf(f, "%f\n", s->pl_cdf + k);
    // if(k >= smooth) { s->pl_cdf[k] /= smooth; for(int i=1;i<smooth;i++) s->pl_cdf[k] += 1.0f/smooth*s->pl_cdf[k-i]; }
    if((s->pl_min + k/(float)s->pl_step < 400) || (s->pl_min + k/(float)s->pl_step > 700)) s->pl_cdf[k] = 0.0f;
  }
  for(int k=1;k<s->pl_num;k++) s->pl_cdf[k] += s->pl_cdf[k-1];
  for(int k=0;k<s->pl_num;k++) s->pl_cdf[k] /= s->pl_cdf[s->pl_num-1];
  s->pl_cdf[s->pl_num-1] = 1.001f;
  fclose(f);

#if 0
  for(int k=0;k<s->pl_num;k+=10) if((s->pl_min + k/(float)s->pl_step >= 400) && (s->pl_min + k/(float)s->pl_step <= 700)) printf("%f %f\n", s->pl_min + k/(float)s->pl_step, s->pl_cdf[k] - s->pl_cdf[k-10]);
  for(int k=1;k<s->pl_num;k++) printf("%d %f %f\n", s->pl_min + s->pl_step * k, s->pl_cdf[k] - s->pl_cdf[k-1], s->pl_cdf[k]);

  // monte carlo experiment:
  float sum = 0.0f;
  for(int i=0;i<10000;i++)
  {
    const float rnd = rand()/(RAND_MAX + 1.0f);
    const float lambda = spectrum_sample_lambda(rnd);
    const float pl = get_spectrum2(s->pl_cdf, s->pl_min, s->pl_step, s->pl_num, lambda);
    sum += pl*(700.0f - 400.0f);
  }
  printf("sum = %f\n", sum);
  sum = 0.0f;
  for(int i=0;i<10000;i++)
  {
    const float rnd = rand()/(RAND_MAX + 1.0f);
    int k = sample_cdf(s->pl_cdf, s->pl_num, rnd);
    const float lambda = s->pl_min + s->pl_step*k;
    sum += 1;//k/(float)s->pl_num;
  }
  printf("sum = %f\n", sum);
#endif

  printf("[medium_fluko] parsed sigma: %d %d, # = %d\n", s->mu_min, s->mu_step, s->mu_num);
  printf("[medium_fluko] parsed pl: %d %d, # = %d, sum = %f\n", s->pl_min, s->pl_step, s->pl_num, s->pl_cdf[s->pl_num-1]);
  return 0;
}
