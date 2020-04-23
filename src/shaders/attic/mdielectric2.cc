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
  float *d = (float *)malloc(2*sizeof(float));
  *data = d;
  int i = fscanf(f, "%f %f", d, d+1);
  if(i < 1)
  {
    fprintf(stderr, "[mdielectric] could not parse all arguments! expecting: n_d [abbe]\n");
    d[0] = 1.5f;
    d[1] = 50;
    return 1;
  }
  if(i != 2) d[1] = 50.0f;
  int dreggn = 0;
  dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) return 1;
  return 0;
}

extern "C" void cleanup(void *data)
{
  free(data);
}

extern "C" float prepare(path_t *p, int v, void *data)
{ 
  float n_d = ((float *)data)[0];
  float V_d = ((float *)data)[1];
  // set volume properties on next segment. we're not using this in sample(), but do it to instruct path space.
  p->v[v].interior.ior = spectrum_eta_from_abbe(n_d, V_d, p->lambda);
  p->v[v].material_modes = static_cast<vertex_scattermode_t>(s_reflect | s_transmit);
  const float eta_ratio = path_eta_ratio(p, v); // n1/n2
  if((p->v[v].shading.roughness > GLOSSY_THR) && (fabsf(eta_ratio - 1.0f) >= 1e-3f)) p->v[v].material_modes = static_cast<vertex_scattermode_t>(p->v[v].material_modes | s_glossy);
  else p->v[v].material_modes = static_cast<vertex_scattermode_t>(p->v[v].material_modes | s_specular);
  return 1.0f;
}

extern "C" float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  float wi[3];
  float wo[3];
  float n[3];
  if(e1 < e2)
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = p->e[e1].omega[k];
      wo[k] = p->e[e2].omega[k];
      n[k]  = p->v[v].hit.n[k];
    }
  }
  else
  {
    for(int k=0;k<3;k++)
    {
      wi[k] = -p->e[e1].omega[k];
      wo[k] = -p->e[e2].omega[k];
      if(p->v[v].mode & s_transmit)
        n[k] = -p->v[v].hit.n[k];
      else
        n[k] = p->v[v].hit.n[k];
    }
  }
  const float wit[3] = {
    dotproduct(p->v[v].hit.a, wi),
    dotproduct(p->v[v].hit.b, wi),
    dotproduct(n, wi)};
  const float wot[3] = {
    dotproduct(p->v[v].hit.a, wo),
    dotproduct(p->v[v].hit.b, wo),
    dotproduct(n, wo)};
  assert(wit[2] <= 0);
  // -wi[2] will always be positive
  //  wo[2] will be positive for R, neg for T

  // if(wit[2] * wot[2] == 0.0) return 0.0f;
  // if(wot[2] > 0.0f && !(p->v[v].mode & s_reflect))  return 0.0f;
  // if(wot[2] < 0.0f && !(p->v[v].mode & s_transmit)) return 0.0f;

  const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(eta_ratio < 0.0f) return 0.0f; // volume nesting broken.
  if(fabsf(eta_ratio - 1.0f) < 1e-3f) // index matched
    return 1.0f;

  MicrosurfaceDielectric surf(true, false, p->v[v].shading.roughness, p->v[v].shading.roughness, 1./eta_ratio);
  return surf.pdf(-vec3(wit[0], wit[1], wit[2]), vec3(wot[0], wot[1], wot[2]));
}


extern "C" float sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.
  const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(eta_ratio < 0.0f) return 0.0f; // volume nesting broken.
  if(fabsf(eta_ratio - 1.0f) < 1e-3f)
  { // index matched
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k];
    p->v[v].mode = static_cast<vertex_scattermode_t>(s_specular | s_transmit);
    p->v[v+1].pdf = 1.0f;
    return p->v[v].shading.rg;
  }

  const float r = p->v[v].shading.roughness;
  if(r > GLOSSY_THR)
  {
    // transform to tangent space
    const float wit[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    MicrosurfaceDielectric surf(true, false, r, r, 1./eta_ratio);
    int bounces = 3;
    vec3 wot;
    float throughput = surf.sample(-vec3(wit[0], wit[1], wit[2]), wot, bounces,
        pointsampler(p, s_dim_scatter_mode),
        pointsampler(p, s_dim_omega_x),
        pointsampler(p, s_dim_omega_y));
    if(throughput <= 0.0) return 0.0f;

    if(wot[2] >= 0) p->v[v].mode = s_reflect;
    else            p->v[v].mode = s_transmit;
    p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy);
    // world space:
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = wot[0]*p->v[v].hit.a[k] + wot[1]*p->v[v].hit.b[k] + wot[2]*p->v[v].hit.n[k];
    p->v[v+1].pdf = surf.pdf(-vec3(wit[0], wit[1], wit[2]), vec3(wot[0], wot[1], wot[2]));
    return p->v[v].shading.rg * throughput;
  }

  // else specular case:
  if(p->v[v].mode & s_reflect)
    p->v[v+1].pdf = .3;// XXX R;
  else
    p->v[v+1].pdf = .7;// XXX 1.0f-R;

  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_specular);
  return p->v[v].shading.rg;
}


extern "C" float brdf(path_t *p, int v, void *data)
{
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(eta_ratio < 0.0f) return 0.0f; // volume nesting broken.
  int index_matched = (fabsf(eta_ratio - 1.0f) < 1e-3f);

  if(cos_out == 0.0f || cos_in == 0.0f) return 0.0f;
  if(!index_matched && (cos_in * cos_out > 0)) p->v[v].mode = s_reflect;
  else p->v[v].mode = s_transmit;

  const float r = p->v[v].shading.roughness;
  if((r > GLOSSY_THR) && !index_matched)
    p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy);
  else
    p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_specular);


  float h[3], cosh = 0.0f;
  if(index_matched)
  {
    const float dot_wo_n = dotproduct(p->e[v+1].omega, p->v[v].hit.n);
    for(int k=0;k<3;k++) h[k] = - p->e[v].omega[k] + p->e[v+1].omega[k] - 2.0f * dot_wo_n * p->v[v].hit.n[k];
    normalise(h);
    cosh = dotproduct(h, p->v[v].hit.n);
    if(cosh < 0.0f) return 0.0f;
    if(cosh < HALFVEC_COS_THR) return 0.0f;
    return p->v[v].shading.rg;
  }
  else if(p->v[v].mode & s_glossy)
  {
    const float wo[3] = {
      dotproduct(p->v[v].hit.a, p->e[v+1].omega),
      dotproduct(p->v[v].hit.b, p->e[v+1].omega),
      dotproduct(p->v[v].hit.n, p->e[v+1].omega)};
    const float wi[3] = {
      dotproduct(p->v[v].hit.a, p->e[v].omega),
      dotproduct(p->v[v].hit.b, p->e[v].omega),
      dotproduct(p->v[v].hit.n, p->e[v].omega)};
    MicrosurfaceDielectric surf(true, false, p->v[v].shading.roughness, p->v[v].shading.roughness, 1./eta_ratio);
    float eval = 0.0f;
    const int num = 1;
    for(int k=0;k<num;k++)
    {
      eval += (fabsf(wi[2]) < fabsf(wo[2])) ?
        surf.eval(-vec3(wi[0], wi[1], wi[2]), vec3(wo[0], wo[1], wo[2]), 0) :
        surf.eval(vec3(wo[0], wo[1], wo[2]), -vec3(wi[0], wi[1], wi[2]), 0);
    }
    return p->v[v].shading.rg * eval/num;
  }

  assert(0); // specular case not implemented
  return 0.0f;
}

