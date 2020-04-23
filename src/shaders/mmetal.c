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
#define MICRO_MATERIAL_CONDUCTOR
#define MICRO_SLOPE_GGX
#include "microfacet.h"
#include "fresnel.h"

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

typedef struct metal_t
{
  int mat;
}
metal_t;

int init(FILE *f, void **data)
{
  metal_t *m = (metal_t *)malloc(sizeof(metal_t));
  char mat[512];
  int i = fscanf(f, "%s", mat);
  if(i < 1)
  {
    fprintf(stderr, "[mmetal] could not parse all arguments! expecting: <ior material name>\n");
    return 1;
  }
  m->mat = fresnel_get_material(mat);
  if(m->mat < 0)
  {
    fprintf(stderr, "[mmetal] WARNING: didn't find `%s' in material list!\n", mat);
    m->mat = 0;
  }
  *data = m;
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) return 1;
  return 0;
}

void cleanup(void *data)
{
  free(data);
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
  const metal_t *m = (metal_t *)data;
  if(!(p->v[v].mode & s_reflect)) return 0.0f;
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

  // XXX
  assert(p->v[v].shading.roughness > 0); // XXX not implemented yet

  const float n1 = p->e[e1].vol.ior;
  const complex float n2 = fresnel_get_ior(m->mat, p->lambda);
  return micro_multiple_pdf(wit, wot,
      p->v[v].shading.roughness, p->v[v].shading.roughness,
      n1, n2, p->index + 1337*v);
}


float sample(path_t *p, void* data)
{
  const metal_t *m = (metal_t *)data;
  const int v = p->length-1; // current vertex.

  const float r = p->v[v].shading.roughness;
  const float n1 = p->e[v].vol.ior;
  const complex float n2 = fresnel_get_ior(m->mat, p->lambda);
  if(r > GLOSSY_THR)
  {
    // wi in tangent space
    float wo[3];
    const float wi[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    float throughput = micro_multiple_sample(wi, wo, r, r, n1, n2,
        p->index + 1337*v,
        pointsampler(p, s_dim_scatter_mode),
        pointsampler(p, s_dim_omega_x),
        pointsampler(p, s_dim_omega_y));

    p->v[v].mode = s_reflect | s_glossy;
    // world space:
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = wo[0]*p->v[v].hit.a[k] + wo[1]*p->v[v].hit.b[k] + wo[2]*p->v[v].hit.n[k];
    p->v[v+1].pdf = micro_multiple_pdf(wi, wo, r, r, n1, n2, p->index + 1337*v);
    return p->v[v].shading.rg * throughput;
  }
  assert(0); // XXX specular not implemented yet
  return 0.0;
}


float brdf(path_t *p, int v, void *data)
{
  const metal_t *m = (metal_t *)data;
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  const float n1 = p->e[v].vol.ior;
  const complex float n2 = fresnel_get_ior(m->mat, p->lambda);

  if(cos_out <= 0.0f || cos_in <= 0.0f) return 0.0f;
  p->v[v].mode = s_reflect | s_glossy; // XXX specular not implemented

  const float r = p->v[v].shading.roughness;
  assert(r > GLOSSY_THR);

  if(p->v[v].mode & s_glossy)
  {
    const float wo[3] = {
      dotproduct(p->v[v].hit.a, p->e[v+1].omega),
      dotproduct(p->v[v].hit.b, p->e[v+1].omega),
      dotproduct(p->v[v].hit.n, p->e[v+1].omega)};
    const float wi[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    return p->v[v].shading.rg * micro_multiple_eval(wi, wo, r, r, n1, n2, p->index + 1337*v);
  }

  assert(0); // specular case not implemented
  return 0.0f;
}

