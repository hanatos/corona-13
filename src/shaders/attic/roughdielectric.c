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
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "shader.h"
#include "spectrum.h"
#include "sampler_common.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define HALFVEC_SQR_HV_EPS 1e-8f
// this is the corresponding cosine threshold: sqrt(1-HALFVEC_SQR_HV_EPS)
// #define HALFVEC_COS_THR .99999999499999998749
// allow some more for weird cases where perfect spheres get voronoi point acne otherwise
#define HALFVEC_COS_THR .999
// phong exponent beyond which specular reflection will be used
#define GLOSSY_THR 10000.0f

int init(FILE *f, void **data)
{
  float *d = (float *)malloc(2*sizeof(float));
  *data = d;
  int i = fscanf(f, "%f %f", d, d+1);
  if(i < 1)
  {
    fprintf(stderr, "[roughdielectric] could not parse all arguments! expecting: n_d [abbe]\n");
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

void cleanup(void *data)
{
  free(data);
}

float prepare(path_t *p, int v, void *data)
{ 
  float n_d = ((float *)data)[0];
  float V_d = ((float *)data)[1];
  // set volume properties on next segment. we're not using this in sample(), but do it to instruct path space.
  p->v[v].interior.ior = spectrum_eta_from_abbe(n_d, V_d, p->lambda);
  p->v[v].material_modes = s_reflect | s_transmit;
  const float roughness = p->v[v].shading.roughness * p->v[v].shading.roughness;
  const float exponent = 2.0f/roughness - 2.0f;
  if((exponent < GLOSSY_THR) && (fabsf(path_eta_ratio(p, v) - 1.0f) >= 1e-3f)) p->v[v].material_modes |= s_glossy;
  else p->v[v].material_modes |= s_specular;
  return 1.0f;
}

static inline float fresnel(const float n1, const float n2, const float cosr, const float cost)
{
  if(cost <= 0.0f) return 1.0f; // total inner reflection
  // fresnel for unpolarized light:
  const float Rs = (n1*cosr - n2*cost)/(n1*cosr + n2*cost);
  const float Rp = (n1*cost - n2*cosr)/(n1*cost + n2*cosr);
  return fminf(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}

static inline float G1(const float *o, const float *n, const float *h, const float roughness)
{
  const float cos_th = fabsf(dotproduct(o, n));
  const float sin_th = sqrtf(1.0f-cos_th*cos_th);
  const float tan_theta = sin_th/cos_th;
  const float a = 1.0f/(roughness * tan_theta);
  if(a < 1.6f) return (3.535*a + 2.181*a*a)/(1.0f + 2.276f*a + 2.577f*a*a);
  else return 1.0f;
}

static inline float shadowing(const float *oi, const float *oo, const float *n, const float *h, const float roughness)
{
  // check sidedness of micro vs macro normal:
  if(dotproduct(oi, n)*dotproduct(oi, h) < 0.0f) return 0.0f;
  if(dotproduct(oo, n)*dotproduct(oo, h) < 0.0f) return 0.0f;
  return G1(oi, n, h, roughness) * G1(oo, n, h, roughness);
}

static inline float roughdie_pdf(path_t *p, int e1, int v, int e2, void *data)
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
  const float cos_in  = -dotproduct(n, wi); // this will always be positive
  const float cos_out =  dotproduct(n, wo); // this will be positive for R, neg for T
  if(cos_in * cos_out == 0.0) return 0.0f;
  if(cos_out > 0.0f && !(p->v[v].mode & s_reflect))  return 0.0f;
  if(cos_out < 0.0f && !(p->v[v].mode & s_transmit)) return 0.0f;
  float h[3];

  // reverse pdf needs different eta ratio
  float n1, n2;
  { // scope braces to make sure eta isn't used after this
    const float eta = path_eta_ratio(p, v); // n1/n2 in tracing direction
    if(eta < 0.0f) return 0.0f; // volume nesting broken.
    if((p->v[v].mode & s_transmit) && (e2 < e1))
    {
      n1 = 1.0;
      n2 = eta;
    }
    else
    {
      n1 = eta;
      n2 = 1.0;
    }
  }

  float cosr = 0.0f, cosh = 0.0f;
  if(fabsf(n1/n2 - 1.0f) < 1e-3f)
  { // index matched
    const float dot_wo_n = dotproduct(wo, n);
    for(int k=0;k<3;k++) h[k] = -wi[k] + wo[k] - 2.0f * dot_wo_n * n[k];
    normalise(h);
    cosh = dotproduct(h, n);
    if(p->v[v].mode != (s_transmit | s_specular)) return 0.0f;
    if(cosh < HALFVEC_COS_THR) return 0.0f;
    return 1.0f;
  }
  else if(p->v[v].mode & s_reflect)
  {
    for(int k=0;k<3;k++) h[k] = wi[k] - wo[k];
    normalise(h);
    cosh = fabsf(dotproduct(h, n));
    cosr = fabsf(dotproduct(h, wi));
  }
  else
  {
    // this h will point into the same hemisphere with the optically thinner medium.
    for(int k=0;k<3;k++) h[k] = n1 * wi[k] - n2 * wo[k];
    normalise(h);
    // we want it to point into the same hemisphere as the tracing direction wi:
    if(n2 < n1) for(int k=0;k<3;k++) h[k] = -h[k];
    cosh = dotproduct(h, n);
    if(cosh < 0.0f) return 0.0f;
    // cosh < 0 means microfacet needs to point down to connect e1 and e2. crazy shit.
    cosr = - dotproduct(h, wi);
    // need to hit micro facet from front:
    if(cosr <= 0.0f) return 0.0f;
  }

  const float cost2 = 1.0f - (n1/n2)*(n1/n2) * (1.0f - cosr*cosr);
  const float cost = (cost2 <= 0.0f) ? 0.0f : sqrtf(cost2);
  const float R = fresnel(n1, n2, cosr, cost);
  const float fudge_factor = 1.2f - 0.2f * sqrtf(fabsf(cos_in));
  const float roughness = p->v[v].shading.roughness * p->v[v].shading.roughness;
  const float exponent_f = 2.0f/(fudge_factor * roughness) - 2.0f;

  float pdf = 1.0f;
  if(p->v[v].mode & s_reflect)
  {
    if(p->v[v].mode & s_specular)//!(exponent < GLOSSY_THR))
    { // pure specular case reflect
      if(cosh < HALFVEC_COS_THR) return 0.0f;
      return R;
    }

    pdf *= 1.0f / (4.0f * fabsf(dotproduct(wo, h)));
    pdf *= R;
  }
  else
  {
    if(p->v[v].mode & s_specular)//!(exponent < GLOSSY_THR))
    { // pure specular case transmit
      if(cosh < HALFVEC_COS_THR) return 0.0f;
      return fmaxf(0.0f, 1.0f-R);
    }
    float denom = n1*cosr - n2*cost;
    pdf *= (n2*n2 * cost) / (denom*denom);
    pdf *= fmaxf(0.0f, 1.-R);
  }

  pdf *= powf(cosh, exponent_f) * (exponent_f + 1.0f)/(2.0f*M_PI);
  pdf /= fabsf(cos_out); // to projected solid angle measure
  // check for degenerate divisions etc:
  if(!(pdf > 0.0f)) return 0.0f;
  return pdf;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  return roughdie_pdf(p, e1, v, e2, data);
}

float sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.
  const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(eta_ratio < 0.0f) return 0.0f; // volume nesting broken.
  if(fabsf(eta_ratio - 1.0f) < 1e-3f)
  { // index matched
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k];
    p->v[v].mode = s_specular | s_transmit;
    p->v[v+1].pdf = 1.0f;
    return p->v[v].shading.rg;
  }
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);

  // from the paper, avoid high weights:
  float fudge_factor = 1.2f - 0.2f * sqrtf(fabsf(cos_in));
  const float roughness = p->v[v].shading.roughness * p->v[v].shading.roughness;
  const float exponent_f = (2.0f/(fudge_factor * roughness) - 2.0f);
  const float exponent   = (2.0f/roughness - 2.0f);

  float *n = p->v[v].hit.n;
  float ht[3] = {0.0, 0.0, 1.0}; // h in tangent space
  float pdf_h = 1.0f;
  if(exponent < GLOSSY_THR)
  {
    sample_cos_k(ht, ht+1, ht+2, exponent_f, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y));
    pdf_h = powf(ht[2], exponent_f) * (exponent_f + 1.0f)/(2.0f*M_PI);
  }
  float pdf = pdf_h;
  float h[3];
  // world space:
  for(int k=0;k<3;k++) h[k] = ht[0]*p->v[v].hit.a[k] + ht[1]*p->v[v].hit.b[k] + ht[2]*n[k];
  const float cosr = -dotproduct(p->e[v].omega, h);
  if(cosr <= 0.0f) return 0.0f; // sampled a microfacet from the wrong side :(

  const float n1 = eta_ratio, n2 = 1.0f; // fake etas, gives same result.

  const float cost2 = 1.0f - eta_ratio*eta_ratio * (1.0f - cosr*cosr);
  const float cost = (cost2 <= 0.0f) ? 0.0f : sqrtf(cost2);
  const float R = fresnel(n1, n2, cosr, cost);

  if(pointsampler(p, s_dim_scatter_mode) < R)
  {
    p->v[v].mode = s_reflect;
    pdf *= R;
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k] + 2.0f * cosr * h[k];
    if(dotproduct(p->e[v+1].omega, n) <= 0.0f) return 0.0f;
    pdf *= 1.0f/(4.0f * cosr);
  }
  else
  { // transmit
    p->v[v].mode = s_transmit;
    pdf *= (1.0-R);
    if(cost2 <= 0.0f) return 0.0f;
    const float f = eta_ratio*cosr - cost;
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k]*eta_ratio + f * h[k];
    normalise(p->e[v+1].omega);
    if(dotproduct(p->e[v+1].omega, n) >= 0.0f) return 0.0f;
    const float denom = n1 * cosr - n2 * cost;
    pdf *= n2 * n2 * cost / (denom*denom);
  }

  if(exponent < GLOSSY_THR)
  {
    p->v[v+1].pdf = pdf/fabsf(dotproduct(p->e[v+1].omega, n)); // convert to projected solid angle measure
    p->v[v].mode |= s_glossy;
  const float roughness = p->v[v].shading.roughness * p->v[v].shading.roughness;
    const float num = shadowing(p->e[v].omega, p->e[v+1].omega, n, h, roughness) * cosr * powf(ht[2], exponent) * (exponent + 1.0f) / (2.0f*M_PI);
    const float den = pdf_h * cos_in;
    if(den == 0.0f) return 0.0f;

    return p->v[v].shading.rg * fabsf(num/den);
  }
  // else specular case:
  if(p->v[v].mode & s_reflect)
    p->v[v+1].pdf = R;
  else
    p->v[v+1].pdf = 1.0f-R;

  p->v[v].mode |= s_specular;
  return p->v[v].shading.rg;
}


float brdf(path_t *p, int v, void *data)
{
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(eta_ratio < 0.0f) return 0.0f; // volume nesting broken.
  const float n1 = eta_ratio, n2 = 1.0f; // fake etas, gives same result.
  int index_matched = (fabsf(eta_ratio - 1.0f) < 1e-3f);

  if(cos_out == 0.0f || cos_in == 0.0f) return 0.0f;
  if(!index_matched && (cos_in * cos_out > 0)) p->v[v].mode = s_reflect;
  else p->v[v].mode = s_transmit;

  const float roughness = p->v[v].shading.roughness * p->v[v].shading.roughness;
  const float exponent = 2.0f/roughness - 2.0f;
  if((exponent < GLOSSY_THR) && !index_matched)
    p->v[v].mode |= s_glossy;
  else
    p->v[v].mode |= s_specular;


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
  else if(p->v[v].mode & s_reflect)
  {
    for(int k=0;k<3;k++) h[k] = - p->e[v].omega[k] + p->e[v+1].omega[k];
    normalise(h);
    cosh = dotproduct(h, p->v[v].hit.n);
    if(cosh < 0.0f) return 0.0f;
  }
  else
  {
    // h will point to thinner medium. if it doesn't, we need to kill the path.
    for(int k=0;k<3;k++) h[k] = n1 * p->e[v].omega[k] - n2 * p->e[v+1].omega[k];
    normalise(h);
    cosh = dotproduct(h, p->v[v].hit.n);
    if(cosh < 0.0f)
    {
      // thinner medium is with wi, but half vector inconsistent:
      if(n1 < n2) return 0.0f;
      for(int k=0;k<3;k++) h[k] = - h[k];
      cosh = - cosh;
    }
    else if (n2 < n1) // thinner medium is n2, but half vector points to n1, too
      return 0.0f;
  }

  // half vector distribution follows phong:
  const float d = powf(cosh, exponent) * (exponent + 1.0f)/(2.0f*M_PI);
  if(d == 0) return 0.0f;

  // fresnel on microfacet:
  const float cosr = -dotproduct(h, p->e[v].omega);
  if(cosr < 0.0f) return 0.0f; // hit microfacet from back side
  const float cost2 = 1.0f - eta_ratio*eta_ratio * (1.0f - cosr*cosr);
  const float cost = (cost2 <= 0.0f) ? 0.0f : sqrtf(cost2);
  // fresnel for unpolarized light:
  const float R = fresnel(n1, n2, cosr, cost);
  const float G = shadowing(p->e[v].omega, p->e[v+1].omega, p->v[v].hit.n, h, roughness);

  if(p->v[v].mode & s_reflect)
  {
    if(cos_in == 0.0f) return 0.0f;
    if(p->v[v].mode & s_glossy)
      return p->v[v].shading.rg * R * d * G / (4.0f*fabsf(cos_in*cos_out));

    // now check dot to outgoing dir for pure specular case
    if(cosh < HALFVEC_COS_THR) return 0.0f;
    if(cost2 < 0.0f) return p->v[v].shading.rg;
    return p->v[v].shading.rg * R;
  }
  else
  {
    const float denom = n1*cosr - n2*cost;
    if(cos_in == 0.0f || denom*denom == 0.0f) return 0.0f;
    if(p->v[v].mode & s_glossy)
      return p->v[v].shading.rg * (1-R) * d * G * n2 * n2 * cosr * cost / (fabsf(cos_in*cos_out) * denom*denom);

    // pure specular case
    if(cosh < HALFVEC_COS_THR) return 0.0f;
    return p->v[v].shading.rg * fmaxf(0.0f, 1.0f-R);
  }
}

