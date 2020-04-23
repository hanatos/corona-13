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
#include <assert.h>

// spectral, reciprocal phong/blinn (diff, spec, mirr, plus fresnel) with measured input data

typedef struct phong_t
{
  int lambda_min, lambda_step, lambda_num;
  float *a, *b, *rd, *rs, *k, *r;
}
phong_t;

float get_spectrum(const phong_t *s, float *v, float lambda)
{
  int l = (lambda - s->lambda_min)/(float)s->lambda_step;
  if(l < 0 || l >= s->lambda_num) return 0.0f;
  return v[l];
}

float powf5(const float x)
{
  const float x2 = x*x;
  return x2*x*x2;
}

float get_fresnel(const float a, const float b, const float cos_theta_in, const float cos_theta_out)
{
  // return fminf(1.0f, a + b*powf5((1.0f - .5f*cos_theta_in)*(1.0f - .5f*cos_theta_out)));
  return fminf(1.0f, a + a*powf5((1.0f - .5f*cos_theta_in)*(1.0f - .5f*cos_theta_out)));
}

float phong_brdf(const float in[3], const float n[3], const float out[3], const float k)
{
  float h[3] = {out[0]-in[0], out[1]-in[1], out[2]-in[2]};
  normalise(h);
  const float dot = fmaxf(0.0f, dotproduct(n, h));
  return (k+1)*powf(dot, k)/(fmaxf(0.001f, -dotproduct(h, in))*fmaxf(-dotproduct(n, in), dotproduct(n, out))*8.0f*M_PI);
}

extern float specularity(const float *omega_in, const rayhit_t *hit, const float rr, void *data)
{
  const phong_t *s = (phong_t *)data;
  const float doti = - dotproduct(hit->normal, omega_in);
  const float f = get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), doti, doti);
  if(rr < f)
  { 
    if(2.0f*rr < f) return 1.0f;
    return 0.5f;//fminf(1.0f, get_spectrum(s, s->k, hit->lambda)/10.0f);
  }
  else return 0.5f;
}

extern float pdf_rr(const float *omega_in, const rayhit_t *hit, const float *omega_out, const float rr, void *data)
{
  const phong_t *s = (phong_t *)data;
  const float doto = dotproduct(hit->normal, omega_out);
  if(doto <= 0) return 0.0f;
  const float doti = - dotproduct(hit->normal, omega_in);

  const float f = get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), doti, doto);
  if(rr < f)
  { 
    if(2.0f*rr < f) return 0.0f;
    float h[3];
    for(int k=0;k<3;k++) h[k] = omega_out[k] - omega_in[k];
    normalise(h);
    const float k = get_spectrum(s, s->k, hit->lambda);
    return (k+1)*powf(dotproduct(hit->normal, h), k)/(8.0f*M_PI*dotproduct(omega_in, h));
  }
  else return doto/M_PI;
}

extern float sample(const float *in, rayhit_t *hit, float *out, const float x1, const float x2, const float rr, void *data)
{
  if(hit->inside) return 0.0f;
  const phong_t *s = (phong_t *)data;

  float cos_theta_in = - dotproduct(hit->normal, in);
  // importance sample by fresnel term (not reciprocal!)
  const float f = get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), cos_theta_in, cos_theta_in);
  // printf("fresnel %f\n", f);
  float x, y, z, h[3];
  if(rr < f)
  {
    // fake reflection
    if(2.0f*rr < f)
    {
      for(int k=0;k<3;k++) out[k] = in[k] + 2*cos_theta_in*hit->normal[k];
      return 2.0f*get_spectrum(s, s->r, hit->lambda);
    }
    // specular part
    const float k = get_spectrum(s, s->k, hit->lambda);
    sample_cos_k(&x, &y, &z, k, x1, x2);
    //printf("k %f\n", k);
    for(int k=0;k<3;k++) h[k] = z*hit->normal[k] + x*hit->a[k] + y*hit->b[k];
    const float hk = - dotproduct(in, h);
    for(int k=0;k<3;k++) out[k] = in[k] + 2*hk*h[k];
    const float cos_theta_out = dotproduct(out, hit->normal);
    if(cos_theta_out <= 0.01f) return 0.0f;
    const float fcorr = get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), cos_theta_in, cos_theta_out)/f;
    return 2.0f*fcorr*get_spectrum(s, s->rs, hit->lambda)/fmaxf(cos_theta_in, cos_theta_out);
  }
  else
  {
    sample_cos(&x, &y, &z, x1, x2);
    for(int k=0; k<3; k++) out[k] = z*hit->normal[k] + x*hit->a[k] + y*hit->b[k];
    const float cos_theta_out = dotproduct(out, hit->normal);
    const float fcorr = (1.0f - get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), cos_theta_in, cos_theta_out))/(1.0f - f);
    return (28.0f/23.0f) * fcorr * get_spectrum(s, s->rd, hit->lambda); // TODO: multiply this at read time.
  }
}

extern float brdf(float *in, rayhit_t *hit, float *out, void *data)
{
  // reflection: p*X = 0.
  const phong_t *s = (phong_t *)data;
  const float cos_theta_in = - dotproduct(in, hit->normal);
  const float cos_theta_out = dotproduct(out, hit->normal);
  if(cos_theta_out <= 0.0f || cos_theta_in <= 0.0f) return 0.0f;
  const float f = get_fresnel(get_spectrum(s, s->a, hit->lambda), get_spectrum(s, s->b, hit->lambda), cos_theta_in, cos_theta_out);
  return (1.0f - f)*28.0*get_spectrum(s, s->rd, hit->lambda)/(23.0f*M_PI) + f*get_spectrum(s, s->rs, hit->lambda)*phong_brdf(in, hit->normal, out, get_spectrum(s, s->k, hit->lambda));
}

extern void cleanup(void *data)
{
  phong_t *s = (phong_t *)data;
  free(s->a); free(s->b); free(s->rd); free(s->rs); free(s->k); free(s->r);
  free(s);
}

extern int init(FILE *f, void **data)
{
  phong_t *s = (phong_t *)malloc(sizeof(phong_t));
  *data = s;
  int num;
  char filename[512];
  if(fscanf(f, "%d %s", &num, filename) != 2)
  {
    printf("[phong] ERROR: could not parse argument! expecting <sample_num> <data.txt>\n");
    return 1;
  }
  FILE *fd = fopen(filename, "rd");
  if(!fd)
  {
    printf("[phong] ERROR: could not open %s!\n", filename);
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(fd, "%d %d %d\n", &(s->lambda_min), &(s->lambda_step), &(s->lambda_num));
  s->a  = (float *)malloc(sizeof(float)*s->lambda_num);
  s->b  = (float *)malloc(sizeof(float)*s->lambda_num);
  s->rd = (float *)malloc(sizeof(float)*s->lambda_num);
  s->rs = (float *)malloc(sizeof(float)*s->lambda_num);
  s->k  = (float *)malloc(sizeof(float)*s->lambda_num);
  s->r  = (float *)malloc(sizeof(float)*s->lambda_num);
  int i[2];
  float v[7];
  int found = 0;
  while (!feof(fd))
  {
    dreggn = fscanf(fd, "%d %d %f %f %f %f %f %f %f", i, i+1, v, v+1, v+2, v+3, v+4, v+5, v+6);
    dreggn = fscanf(fd, "%*[^\n]\n");
    if(i[0] == num)
    {
      found = 1;
      // error [i[1]] = v[0];
      s->a [i[1]] = v[1];
      s->b [i[1]] = v[2];
      s->rd[i[1]] = v[3];
      s->rs[i[1]] = v[4];
      s->k [i[1]] = v[5];
      s->r [i[1]] = v[6];
      // printf("a %f b %f rd %f rs %f k %f r %f\n", s->a[i[1]], s->b[i[1]], s->rd[i[1]], s->rs[i[1]], s->k[i[1]], s->r[i[1]]);
    }
  }
  fclose(fd);
  if(!found)
  {
    fprintf(stderr, "[phong] ERROR: material %d not in file %s!\n", num, filename);
    cleanup(s);
    return 1;
  }

  dreggn = fscanf(f, "%*[^\n]\n");
  return 0;
}


