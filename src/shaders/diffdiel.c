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
#include "sampler_common.h"
#include "shader.h"
#include "spectrum.h"
#include "ggx.h"
#include "points.h"

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

int init(FILE *f, void **data)
{
  float *d = malloc(2*sizeof(float));
  *data = d;
  int i = fscanf(f, "%f %f", d, d+1);
  if(i < 1)
  {
    fprintf(stderr, "[diffdiel] could not parse all arguments! expecting: n_d [abbe]\n");
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

static inline int indexmatched(mf_t n1, mf_t n2)
{
  const mf_t eta = mf_div(n1, n2);
  return mf_any(mf_lt(mf_abs(mf_sub(mf_set1(1.0f), eta)), mf_set1(1e-3f)));
}

float prepare(path_t *p, int v, void *data)
{ 
  float n_d = ((float *)data)[0];
  float V_d = ((float *)data)[1];
  // set volume properties on next segment. we're not using this in sample, but do it to instruct path space.
  p->v[v].interior.ior = spectrum_eta_from_abbe(n_d, V_d, p->lambda);
  p->v[v].material_modes = s_reflect | s_transmit;

  mf_t eta = path_eta_ratio(p, v); // need symmetric degeneration:
  if(indexmatched(eta, mf_set1(1.0f))) p->v[v].shading.roughness = 0.0f;
  if(p->v[v].shading.roughness > GLOSSY_THR) p->v[v].material_modes |= s_glossy;
  // if((p->v[v].shading.roughness > GLOSSY_THR) && (fabsf(path_eta_ratio(p, v) - 1.0f) >= 1e-3f)) p->v[v].material_modes |= s_glossy;
  else p->v[v].material_modes |= s_specular;
  return 1.0f;
}

static inline mf_t fresnel(const mf_t n1, const mf_t n2, const mf_t cosr, const mf_t cost)
{
  // fresnel for unpolarized light:
  const mf_t r1 = mf_mul(n1, cosr), r2 = mf_mul(n2, cosr),
             t1 = mf_mul(n1, cost), t2 = mf_mul(n2, cost);
  const mf_t Rs = mf_div(mf_sub(r1, t2), mf_add(r1, t2));
  const mf_t Rp = mf_div(mf_sub(t1, r2), mf_add(t1, r2));
  return mf_select(
      mf_set1(1.0f), // return R=1 for total internal reflection
      mf_clamp((Rs*Rs + Rp*Rp)*.5f, 0.0f, 1.0f),
      mf_lte(cost, mf_set1(0.0f)));
}

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
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
  if(cos_in * cos_out == 0.0f) return mf_set1(0.0f);
  if(cos_out > 0.0f && !(p->v[v].mode & s_reflect))  return mf_set1(0.0f);
  if(cos_out < 0.0f && !(p->v[v].mode & s_transmit)) return mf_set1(0.0f);
  float h[3];

  // reverse pdf needs different eta ratio
  mf_t n1, n2;
  { // scope braces to make sure eta isn't used after this
    const mf_t eta = path_eta_ratio(p, v); // n1/n2 in tracing direction
    if(mf_all(mf_lt(eta,  mf_set1(0.0f)))) return mf_set1(0.0f); // volume nesting broken.
    if((p->v[v].mode & s_transmit) && (e2 < e1))
    {
      n1 = mf_set1(1.0f);
      n2 = eta;
    }
    else
    {
      n1 = eta;
      n2 = mf_set1(1.0f);
    }
  }

  // need mask instead of pretty much all return statements here:
  mf_t mask = mf_set1(0.0f); // start with nothing masked out

  mf_t cosr = mf_set1(0.0f), cosh = mf_set1(0.0f);
  if(indexmatched(n1, n2))
  { // index matched
    const float dot_wo_n = dotproduct(wo, n);
    for(int k=0;k<3;k++) h[k] = -wi[k] + wo[k] - 2.0f * dot_wo_n * n[k];
    normalise(h);
    cosh = mf_set1(dotproduct(h, n));
    if(p->v[v].mode != (s_transmit | s_specular)) return mf_set1(0.0f);
    if(mf(cosh, 0) < HALFVEC_COS_THR) return mf_set1(0.0f);
    return mf_set1(1.0f);
  }
  else if(p->v[v].mode & s_reflect)
  {
    for(int k=0;k<3;k++) h[k] = wi[k] - wo[k];
    normalise(h);
    cosh = mf_set1(fabsf(dotproduct(h, n)));
    cosr = mf_set1(fabsf(dotproduct(h, wi)));
  }

  // non reciprocal fresnel on cos_in
  const mf_t nr = mf_div(n1, n2);
  const mf_t cost2 = mf_sub(mf_set1(1.0f), mf_mul(mf_mul(nr, nr), mf_sub(mf_set1(1.0f), mf_set1(cos_in*cos_in))));
  const mf_t cost = mf_select(mf_set1(0.0f), mf_sqrt(cost2), mf_lte(cost2, mf_set1(0.0f)));
  mf_t R = fresnel(n1, n2, mf_set1(cos_in), cost);
  // do the precise same culling dance as for sampling here
  if(p->v[v].culled_modes & s_reflect)
  {
    mask = mf_or(mask, mf_eq(R, mf_set1(1.0f))); // mask out everything with total internal reflection
    R = mf_set1(0.0f);
  }
  if(p->v[v].culled_modes & s_transmit)
  {
    mask = mf_or(mask, mf_eq(R, mf_set1(0.0f))); // mask out everything without reflection
    R = mf_set1(1.0f);
  }

  mf_t pdf = mf_set1(1.0f);
  if(p->v[v].mode & s_reflect)
  {
    if(p->v[v].mode & s_specular)
    { // pure specular case reflect
      mask = mf_or(mask, mf_lt(cosh, mf_set1(HALFVEC_COS_THR)));
      return mf_select(mf_set1(0.0f), R, mask);
    }

    pdf = mf_mul(pdf, mf_set1(1.0f / (4.0f * fabsf(dotproduct(wo, h)))));
    pdf = mf_mul(pdf, R);
  }
  else
  {
    if(p->v[v].mode & s_specular)
    { // pure specular case transmit
      mask = mf_or(mask, mf_lt(cosh, mf_set1(HALFVEC_COS_THR)));
      return mf_select(mf_set1(0.0f),
          mf_clamp(mf_sub(mf_set1(1.0f), R), 0.0f, 1.0f),
          mask);
    }
    pdf = mf_set1(1.0f/M_PI);
    pdf = mf_mul(pdf, mf_clamp(mf_sub(mf_set1(1.0f), R), 0.0f, 1.0f));
    mask = mf_or(mask, mf_not(mf_gt(pdf, mf_set1(0.0f))));
    return mf_select(mf_set1(0.0f), pdf, mask);
  }

  pdf = mf_mul(pdf, ggx_pdf_h_mf(cosh, mf_set1(cos_in), cosr, p->v[v].shading.roughness));
  pdf = mf_div(pdf, mf_set1(fabsf(cos_out))); // to projected solid angle measure
  // mask out for degenerate divisions etc:
  mask = mf_or(mask, mf_not(mf_gt(pdf, mf_set1(0.0f))));
  return mf_select(mf_set1(0.0f), pdf, mask);
}


mf_t sample(path_t *p, void* data)
{
  const int v = p->length-1; // current vertex.
  const mf_t eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(mf(eta_ratio, 0) < 0.0f) return mf_set1(0.0f); // volume nesting broken.
  if(indexmatched(eta_ratio, mf_set1(1.0f)))
  { // index matched
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k];
    p->v[v].mode = s_specular | s_transmit;
    p->v[v+1].pdf = mf_set1(1.0f);
    return p->v[v].shading.rg;
  }

  float *n = p->v[v].hit.n;
  float ht[3] = {0.0, 0.0, 1.0}; // h in tangent space
  float pdf_h = 1.0f;
  float h[3] = {n[0], n[1], n[2]}; // init for smooth dielectric
  const float r = p->v[v].shading.roughness;
  const float cos_in = - dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float rand1 = pointsampler(p, s_dim_omega_x);
  const float rand2 = pointsampler(p, s_dim_omega_y);

  if(r > GLOSSY_THR)
  {
    // flip incoming direction to point away from intersection point
    const float wit[3] = {
      - dotproduct(p->v[v].hit.a, p->e[v].omega),
      - dotproduct(p->v[v].hit.b, p->e[v].omega),
      cos_in};
    ggx_sample_h(wit, r, r, rand1, rand2, ht);
    // world space:
    for(int k=0;k<3;k++) h[k] = ht[0]*p->v[v].hit.a[k] + ht[1]*p->v[v].hit.b[k] + ht[2]*n[k];
    pdf_h = ggx_pdf_h(p->e[v].omega, h, n, r);
  }
  float pdf = pdf_h;

  const float cosr = -dotproduct(p->e[v].omega, h);
  if(cosr <= 0.0f) return mf_set1(0.0f); // sampled a microfacet from the wrong side :(

  // non-reciprocal fake fresnel:
  const mf_t n1 = eta_ratio, n2 = mf_set1(1.0f); // fake etas, gives same result.
  const mf_t nr = mf_div(n1, n2);
  const mf_t cost2 = mf_sub(mf_set1(1.0f), mf_mul(mf_mul(nr, nr), mf_sub(mf_set1(1.0f), mf_set1(cos_in*cos_in))));
  const mf_t cost = mf_select(mf_set1(0.0f), mf_sqrt(cost2), mf_lte(cost2, mf_set1(0.0f)));
  mf_t R = fresnel(n1, n2, mf_set1(cos_in), cost);
  mf_t mask = mf_set1(0.0f);
  if(p->v[v].culled_modes & s_reflect)
  {
    mask = mf_or(mask, mf_eq(R, mf_set1(1.0f))); // mask out everything with total internal reflection
    R = mf_set1(0.0f);
  }
  if(p->v[v].culled_modes & s_transmit)
  {
    mask = mf_or(mask, mf_eq(R, mf_set1(0.0f))); // mask out everything without reflection
    R = mf_set1(1.0f);
  }

  if(pointsampler(p, s_dim_scatter_mode) <= mf(R, 0))
  {
    p->v[v].mode = s_reflect;
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k] + 2.0f * cosr * h[k];
    if(dotproduct(p->e[v+1].omega, n) <= 0.0f) return mf_set1(0.0f);
    pdf *= 1.0f/(4.0f * cosr);

    if(r > GLOSSY_THR)
    {
      p->v[v+1].pdf = mf_mul(R, mf_set1(pdf/fabsf(dotproduct(p->e[v+1].omega, n)))); // convert to projected solid angle measure
      p->v[v].mode |= s_glossy;
      if(dotproduct(p->e[v+1].omega, n)*dotproduct(p->e[v+1].omega, h) < 0.0f) return mf_set1(0.0f);
      return mf_select(mf_set1(0.0f), 
          mf_mul(p->v[v].shading.rg, mf_set1(ggx_shadowing_smith_G1(p->e[v+1].omega, n, h, p->v[v].shading.roughness))),
          mask);
    }

    // else: specular:
    p->v[v+1].pdf = mf_select(mf_set1(0.0f), R, mask);
    p->v[v].mode |= s_specular;
    return mf_select(mf_set1(0.0f), p->v[v].shading.rg, mask);
  }
  else
  { // transmit
    p->v[v].mode = s_glossy | s_transmit;
    p->v[v+1].pdf = mf_div(mf_sub(mf_set1(1.0f), R), mf_set1(M_PI));
    float wo[3];
    sample_cos(wo, wo+1, wo+2, rand1, rand2);
    for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->v[v].hit.a[k]*wo[0] + p->v[v].hit.b[k]*wo[1] - p->v[v].hit.n[k]*wo[2];
    return p->v[v].shading.rg;
  }
}


mf_t brdf(path_t *p, int v, void *data)
{
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  const mf_t eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(mf(eta_ratio, 0) < 0.0f) return mf_set1(0.0f); // volume nesting broken.
  const mf_t n1 = eta_ratio, n2 = mf_set1(1.0f); // fake etas, gives same result.
  int index_matched = indexmatched(n1, n2);

  if(cos_out == 0.0f || cos_in == 0.0f) return mf_set1(0.0f);
  if(!index_matched && (cos_in * cos_out > 0)) p->v[v].mode = s_reflect;
  else p->v[v].mode = s_transmit;

  const float r = p->v[v].shading.roughness;
  if((r > GLOSSY_THR) && !index_matched)
    p->v[v].mode |= s_glossy;
  else
    p->v[v].mode |= s_specular;

  // non-reciprocal fresnel on surface normal
  const mf_t nr = mf_div(n1, n2);
  const mf_t cost2 = mf_sub(mf_set1(1.0f), mf_mul(mf_mul(nr, nr), mf_sub(mf_set1(1.0f), mf_mul(mf_set1(cos_in), mf_set1(cos_in)))));
  const mf_t cost = mf_select(mf_set1(0.0f), mf_sqrt(cost2), mf_lte(cost2, mf_set1(0.0f)));
  mf_t R = fresnel(n1, n2, mf_set1(cos_in), cost);

  float h[3], cosh = 0.0f;
  if(index_matched)
  {
    const float dot_wo_n = dotproduct(p->e[v+1].omega, p->v[v].hit.n);
    for(int k=0;k<3;k++) h[k] = - p->e[v].omega[k] + p->e[v+1].omega[k] - 2.0f * dot_wo_n * p->v[v].hit.n[k];
    normalise(h);
    cosh = dotproduct(h, p->v[v].hit.n);
    if(cosh < 0.0f) return mf_set1(0.0f);
    if(cosh < HALFVEC_COS_THR) return mf_set1(0.0f);
    return p->v[v].shading.rg;
  }
  else if(p->v[v].mode & s_reflect)
  {
    for(int k=0;k<3;k++) h[k] = - p->e[v].omega[k] + p->e[v+1].omega[k];
    normalise(h);
    cosh = dotproduct(h, p->v[v].hit.n);
    if(cosh < 0.0f) return mf_set1(0.0f);
    float cosr = dotproduct(h, p->e[v+1].omega);

    // ggx half vector distribution:
    const float DG1 = (p->v[v].mode & s_specular) ? 1.0f : ggx_pdf_h(p->e[v].omega, h, p->v[v].hit.n, p->v[v].shading.roughness);
    if(DG1 == 0) return mf_set1(0.0f);
    // check microfacet/macrofacet orientation consistency (does not make any difference, seems to be good already)
    // if(dotproduct(p->e[v].omega, p->v[v].hit.n)*dotproduct(p->e[v].omega, h) < 0.0f) return 0.0f;
    // if(dotproduct(p->e[v+1].omega, p->v[v].hit.n)*dotproduct(p->e[v+1].omega, h) < 0.0f) return 0.0f;

    const float G1 = ggx_shadowing_smith_G1(p->e[v+1].omega, p->v[v].hit.n, h, p->v[v].shading.roughness);
    if(cos_in == 0.0f) return mf_set1(0.0f);
    if(p->v[v].mode & s_glossy)
      return mf_mul(mf_mul(p->v[v].shading.rg, R), mf_set1(DG1 * G1 / (4.0f*fabsf(cosr*cos_out))));

    // now check dot to outgoing dir for pure specular case
    if(cosh < HALFVEC_COS_THR) return mf_set1(0.0f);
    return mf_mul(p->v[v].shading.rg, R);
  }
  else
  {
    if(cos_in == 0.0f) return mf_set1(0.0f);
    if(p->v[v].mode & s_glossy)
      return mf_mul(p->v[v].shading.rg, mf_div(mf_clamp(mf_sub(mf_set1(1.0f), R), 0.0f, 1.0f), mf_set1(M_PI)));

    // pure specular case
    if(cosh < HALFVEC_COS_THR) return mf_set1(0.0f);
    return mf_mul(p->v[v].shading.rg, mf_clamp(mf_sub(mf_set1(1.0f), R), 0.0f, 1.0f));
  }
}

