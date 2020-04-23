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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float sample(const float* omega_in, rayhit_t* hit, float* omega_out, const float qx, const float qy, const float rr, void* data)
{
  const float s = *(float*)&data;
  const float dot = hit->normal[0]*omega_in[0] + hit->normal[1]*omega_in[1] + hit->normal[2]*omega_in[2];
  const float fresnel = s + (1.0f - s)*powf((1.0f - fabsf(dot)), 5.0f);
  const float p_s = fresnel;//*hit->rs;
  // const float p_d = (1.0f - p_s)*hit->rd;
  if(rr < p_s)
  {
    for(int k=0;k<3;k++) omega_out[k] = omega_in[k] - 2*dot*hit->normal[k];
    return hit->rs;//1.0f;
  }
  else// if(rr < p_s + p_d)
  {
    for(int k=0;k<3;k++) omega_out[k] = omega_in[k];
    return hit->rd;//1.0f;
  }
  //else return 0.0f;
}

float pdf(const float *omega_in, const rayhit_t *hit, const float *omega_out, void *data)
{
  return -1.0f;
}

float brdf(path_t *p, int v, void *data)
{
  return 0.0f;
}

void cleanup() {}

int init(FILE* f, void** data)
{
  float *d = (float*)data;
  int i = fscanf(f, "%f", d);
  if(i != 1)
  {
    fprintf(stderr, "[glass] could not parse all arguments! expecting glass specular\n");
    *d = 0.2f;
    return 1;
  }
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  return 0;
}
