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

extern "C" {
#include "corona_common.h"
#include "shader.h"
#include "spectrum_common.h"
}
#include "MicrosurfaceScattering.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define HALFVEC_SQR_HV_EPS 1e-8f
// this is the corresponding cosine threshold: sqrt(1-HALFVEC_SQR_HV_EPS)
// #define HALFVEC_COS_THR .99999999499999998749
// allow some more for weird cases where perfect spheres get voronoi point acne otherwise
#define HALFVEC_COS_THR .999
// roughness below which specular reflection will be used
#define GLOSSY_THR 1e-4f

extern "C" int init(FILE *f, void **data)
{
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) return 1;
  return 0;
}

extern "C" float prepare(path_t *p, int v, void *data)
{ 
  p->v[v].material_modes = static_cast<vertex_scattermode_t>(s_reflect | s_transmit);
  p->v[v].material_modes = static_cast<vertex_scattermode_t>(p->v[v].material_modes | s_glossy);
  return 1.0f;
}

extern "C" float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  // let's lie a little here:
  return 1.0/M_PI;
}


extern "C" float sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.
  const float r = p->v[v].shading.roughness;
  // wi in tangent space
  const float wit[3] = {
    dotproduct(p->v[v].hit.a, p->e[v].omega),
    dotproduct(p->v[v].hit.b, p->e[v].omega),
    dotproduct(p->v[v].hit.n, p->e[v].omega)};

  MicrosurfaceDiffuse surf(true, false, r, r, p->v[v].shading.rd);
  int bounces = 3;
  vec3 wot;
  float throughput = surf.sample(-vec3(wit[0], wit[1], wit[2]), wot, bounces,
      pointsampler(p, s_dim_scatter_mode),
      pointsampler(p, s_dim_omega_x),
      pointsampler(p, s_dim_omega_y));
  if(throughput <= 0.0) return 0.0f;

  p->v[v].mode = static_cast<vertex_scattermode_t>(s_reflect | s_glossy);
  // world space:
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = wot[0]*p->v[v].hit.a[k] + wot[1]*p->v[v].hit.b[k] + wot[2]*p->v[v].hit.n[k];
  p->v[v+1].pdf = 1.0/M_PI;
  return throughput;
}


extern "C" float brdf(path_t *p, int v, void *data)
{
  const float r = p->v[v].shading.roughness;
  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy);

  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  if(cos_out <= 0.0f || cos_in <= 0.0f) return 0.0f;

  const float wo[3] = {
    dotproduct(p->v[v].hit.a, p->e[v+1].omega),
    dotproduct(p->v[v].hit.b, p->e[v+1].omega),
    dotproduct(p->v[v].hit.n, p->e[v+1].omega)};
  const float wi[3] = {
    dotproduct(p->v[v].hit.a, p->e[v].omega),
    dotproduct(p->v[v].hit.b, p->e[v].omega),
    dotproduct(p->v[v].hit.n, p->e[v].omega)};
  MicrosurfaceDiffuse surf(true, false, r, r, p->v[v].shading.rd);
  float eval = 0.0f;
  const int num = 1;
  for(int k=0;k<num;k++)
    eval += (fabsf(wi[2]) < fabsf(wo[2])) ?
      surf.eval(-vec3(wi[0], wi[1], wi[2]), vec3(wo[0], wo[1], wo[2]), 0) :
      surf.eval(vec3(wo[0], wo[1], wo[2]), -vec3(wi[0], wi[1], wi[2]), 0);
  return eval/num;
}

