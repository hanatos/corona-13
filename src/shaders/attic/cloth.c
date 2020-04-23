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

float sample(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
#if 0
  float x, y, z;
  sample_cos(&x, &y, &z, qx, qy);
  for(int k=0;k<3;k++) omega_out[k] = hit->a[k]*x + hit->b[k]*y;
  if(rr < 0.5f) for(int k=0;k<3;k++) omega_out[k] += hit->normal[k]*z;
  else for(int k=0;k<3;k++) omega_out[k] -= hit->normal[k]*z;
  return hit->rd;
#else
  float s = *(float*)&data;
  const float dot = dotproduct(hit->normal, omega_in);
  const float fresnel = hit->rs + (1.0f - hit->rs)*powf((1.0f - fabsf(dot)), 5.0f);
  float x, y, z;
  sample_cos_k(&x, &y, &z, s, qx, qy);
  float a[3], b[3], r[3];
  if(rr < fresnel) for(int k=0;k<3;k++) r[k] = omega_in[k] - 2.0f * dot * hit->normal[k];
  else for(int k=0;k<3;k++) r[k] = omega_in[k];
  get_onb(r, a, b);
  for(int k=0;k<3;k++) omega_out[k] = z*r[k] + x*a[k] + y*b[k];
  // return (2.0f*M_PI)/(s + 1.0f) * hit->rd;
  return hit->rd;
#endif
}

float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return 1.0f/M_PI;
}

float brdf(float *omega_in, rayhit_t *hit, float *omega_out, void *data)
{
#if 0
  const float dot = dotproduct(hit->normal, omega_out);
  return hit->rd*fabsf(dot)/M_PI;
#else
  const float exp = *(float*)&data;
  float r1[3];
  const float dot = dotproduct(hit->normal, omega_in);
  const float fresnel = hit->rs + (1.0f - hit->rs)*powf((1.0f - fabsf(dot)), 5.0f);
  for(int k=0;k<3;k++) r1[k] = omega_in[k] - 2.0f * dot * hit->normal[k];
  const double dot1 = fmaxf(0.0f, dotproduct(r1, omega_out));
  const double dot2 = fmaxf(0.0f, dotproduct(omega_in, omega_out));
  return hit->rd*(fresnel*((exp + 1.0)*powf(dot1, exp))/(2.0*M_PI) + 
         (1.0f - fresnel)*((exp + 1.0)*powf(dot2, exp))/(2.0*M_PI));
#endif
}

// avoid free(data)
extern void cleanup() {}

extern int init(FILE* f, void** data)
{
  float *d = (float*)data;
  int i = fscanf(f, "%f", d);
  if(i != 1)
  {
    fprintf(stderr, "[cloth::init]: could not parse all arguments! expecting specular\n");
    *d = 0.2f;
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  return 0;
}
