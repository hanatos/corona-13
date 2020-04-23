/*
    This file is part of corona-13.
    copyright (c) 2015 johannes hanika.

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
#define MICRO_MATERIAL_DIFFUSE
#define MICRO_SLOPE_GGX
#include "microfacet.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

// roughness below which specular reflection will be used
#define GLOSSY_THR 1e-4f

int init(FILE *f, void **data)
{
  int dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) return 1;
  return 0;
}

float prepare(path_t *p, int v, void *data)
{ 
  p->v[v].material_modes = s_reflect;
  if(p->v[v].shading.roughness > GLOSSY_THR) p->v[v].material_modes |= s_glossy;
  else p->v[v].material_modes |= s_specular;
  return 1.0f;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_reflect)) return 0.0f;
  if(p->v[v].shading.roughness == 0.0) return 1.0f/M_PI; // smooth diffuse
  float wi[3];
  float wo[3];
  if(e1 < e2)
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = p->e[e1].omega[k];
      wo[k] = p->e[e2].omega[k];
    }
  }
  else
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = -p->e[e1].omega[k];
      wo[k] = -p->e[e2].omega[k];
    }
  }
  const float wit[3] = {
    dotproduct(p->v[v].hit.a, wi),
    dotproduct(p->v[v].hit.b, wi),
    dotproduct(p->v[v].hit.n, wi)};
  const float wot[3] = {
    dotproduct(p->v[v].hit.a, wo),
    dotproduct(p->v[v].hit.b, wo),
    dotproduct(p->v[v].hit.n, wo)};

  if(wit[2] >= 0.0f) return 0.0f;
  if(wot[2] <= 0.0f) return 0.0f;

  return micro_multiple_pdf(wit, wot,
      p->v[v].shading.roughness, p->v[v].shading.roughness,
      p->v[v].shading.rd, 0.0, p->index + 1337*v);
}


float sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.

  const float r = p->v[v].shading.roughness;
  float wo[3];
  if(r > GLOSSY_THR)
  { // rough diffuse
    // wi in tangent space
    const float wi[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    float throughput = micro_multiple_sample(wi, wo, r, r, p->v[v].shading.rd, 0.0,
        p->index + 1337*v,
        pointsampler(p, s_dim_scatter_mode),
        pointsampler(p, s_dim_omega_x),
        pointsampler(p, s_dim_omega_y));

    p->v[v].mode = s_reflect | s_glossy;
    // world space:
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = wo[0]*p->v[v].hit.a[k] + wo[1]*p->v[v].hit.b[k] + wo[2]*p->v[v].hit.n[k];
    p->v[v+1].pdf = micro_multiple_pdf(wi, wo, r, r, p->v[v].shading.rd, 0.0, p->index + 1337*v);
    return throughput;
  }
  else
  { // smooth diffuse (no shading normal madness implemented here)
    sample_cos(wo, wo+1, wo+2, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y));
    p->v[v].mode = s_reflect | s_diffuse;
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = wo[0]*p->v[v].hit.a[k] + wo[1]*p->v[v].hit.b[k] + wo[2]*p->v[v].hit.n[k];
    p->v[v+1].pdf = 1.0f/M_PI;
    return p->v[v].shading.rd;
  }
  return 0.0;
}


float brdf(path_t *p, int v, void *data)
{
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  if(cos_out <= 0.0f || cos_in <= 0.0f) return 0.0f;
  p->v[v].mode = s_reflect;
  const float r = p->v[v].shading.roughness;

  if(r > GLOSSY_THR)
  {
    const float wo[3] = {
      dotproduct(p->v[v].hit.a, p->e[v+1].omega),
      dotproduct(p->v[v].hit.b, p->e[v+1].omega),
      dotproduct(p->v[v].hit.n, p->e[v+1].omega)};
    const float wi[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    p->v[v].mode |= s_glossy;
    return micro_multiple_eval(wi, wo, r, r, p->v[v].shading.rd, 0.0, p->index + 1337*v);
  }
  // smooth diffuse
  p->v[v].mode |= s_diffuse;
  return p->v[v].shading.rd / M_PI;
}

