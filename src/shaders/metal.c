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
#include "ggx.h"
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
    fprintf(stderr, "[metal] could not parse all arguments! expecting: <ior material name>\n");
    return 1;
  }
  m->mat = fresnel_get_material(mat);
  if(m->mat < 0)
  {
    fprintf(stderr, "[metal] WARNING: didn't find `%s' in material list!\n", mat);
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

static inline mf_t fresnel(const mf_t n1, const mf_t n2, const mf_t k2, const float cosr)
{
  // compute eta^2 from n1 and complex n2:
  mf_t etar =        mf_div(mf_mul(n1, n2), mf_add(mf_mul(n2, n2), mf_mul(k2, k2)));
  mf_t etai = mf_neg(mf_div(mf_mul(n1, k2), mf_add(mf_mul(n2, n2), mf_mul(k2, k2))));
  mf_t eta2r = mf_sub(mf_mul(etar, etar), mf_mul(etai, etai));
  mf_t eta2i = mf_mul(mf_mul(mf_set1(2.0f), etar), etai);

  // float c_n1 = mf(n1,0);
  // complex float c_n2 = mf(n2,0) + I*mf(k2,0);
  // complex float c_eta2 = c_n1/c_n2*c_n1/c_n2;
  // assert(mf(eta2r,0) - crealf(c_eta2) < 0.001f);
  // assert(mf(eta2i,0) - cimagf(c_eta2) < 0.001f);

  const float sinr = 1.0f - cosr*cosr;
  // compute transmitted cosine cost (complex valued):
  mf_t cost2r = mf_sub(mf_set1(1.0f), mf_mul(eta2r, mf_set1(sinr)));
  mf_t cost2i = mf_mul(eta2i, mf_set1(- sinr));
  mf_t len = mf_sqrt(mf_add(mf_mul(cost2r, cost2r), mf_mul(cost2i, cost2i)));
  mf_t costr = mf_sqrt(mf_mul(mf_set1(0.5f), mf_add(cost2r, len)));
  mf_t costi = mf_sqrt(mf_mul(mf_set1(0.5f), mf_sub(len, cost2r)));
  costi = mf_select(mf_neg(costi), costi, mf_lt(cost2i, mf_set1(0.0f)));

  // const complex float c_cost = csqrtf(1.0f - c_eta2 * (1.0f - cosr*cosr));
  // assert(fabsf(crealf(c_cost) - mf(costr, 0)) < 0.001f);
  // assert(fabsf(cimagf(c_cost) - mf(costi, 0)) < 0.001f);
  
  // compute all mixed terms of n[12] * cos[rt]
  mf_t n1cosr  = mf_mul(n1, mf_set1(cosr));
  mf_t n2cosrr = mf_mul(n2, mf_set1(cosr));
  mf_t n2cosri = mf_mul(k2, mf_set1(cosr));
  mf_t n1costr = mf_mul(n1, costr);
  mf_t n1costi = mf_mul(n1, costi);
  mf_t n2costr = mf_sub(mf_mul(n2, costr), mf_mul(k2, costi));
  mf_t n2costi = mf_add(mf_mul(k2, costr), mf_mul(n2, costi));
  // since we're only interested in the absolute value, we can do
  // the division on the absolute values of the operands:

  // assert(fabsf(crealf(c_n1*cosr) - mf(n1cosr,0)) < 0.001f);
  // assert(fabsf(crealf(c_n2*cosr) - mf(n2cosrr,0)) < 0.001f);
  // assert(fabsf(cimagf(c_n2*cosr) - mf(n2cosri,0)) < 0.001f);
  // assert(fabsf(crealf(c_n1*c_cost) - mf(n1costr,0)) < 0.001f);
  // assert(fabsf(cimagf(c_n1*c_cost) - mf(n1costi,0)) < 0.001f);
  // assert(fabsf(crealf(c_n2*c_cost) - mf(n2costr,0)) < 0.001f);
  // assert(fabsf(cimagf(c_n2*c_cost) - mf(n2costi,0)) < 0.001f);


  mf_t Rs_num = (mf_add(
        mf_mul(mf_sub(n1cosr, n2costr), mf_sub(n1cosr, n2costr)),
        mf_mul(n2costi, n2costi)
        ));
  mf_t Rs_den = (mf_add(
        mf_mul(mf_add(n1cosr, n2costr), mf_add(n1cosr, n2costr)),
        mf_mul(n2costi, n2costi)
        ));
  mf_t Rs2 = mf_div(Rs_num, Rs_den);

  // const float c_Rs_num = cabsf(c_n1*  cosr - c_n2*c_cost);
  // const float c_Rs_den = cabsf(c_n1*  cosr + c_n2*c_cost);
  // const float c_Rs = cabsf((c_n1*  cosr - c_n2*c_cost)/(c_n1*  cosr + c_n2*c_cost));
  // assert(fabsf(c_Rs_num*c_Rs_num - mf(Rs_num,0)) < 0.005f);
  // assert(fabsf(c_Rs_den*c_Rs_den - mf(Rs_den,0)) < 0.005f);
  // assert(fabsf(c_Rs*c_Rs - mf(Rs2,0)) < 0.005f);

  mf_t Rp_num = (mf_add(
        mf_mul(mf_sub(n1costr, n2cosrr), mf_sub(n1costr, n2cosrr)),
        mf_mul(mf_sub(n1costi, n2cosri), mf_sub(n1costi, n2cosri))
        ));
  mf_t Rp_den = (mf_add(
        mf_mul(mf_add(n1costr, n2cosrr), mf_add(n1costr, n2cosrr)),
        mf_mul(mf_add(n1costi, n2cosri), mf_add(n1costi, n2cosri))
        ));
  mf_t Rp2 = mf_div(Rp_num, Rp_den);

  // const float c_Rp = cabsf((c_n1*c_cost - c_n2*  cosr)/(c_n1*c_cost + c_n2*  cosr));
  // assert(fabsf(c_Rp*c_Rp - mf(Rp2,0)) < 0.004f);
      
  return mf_clamp(mf_mul(mf_add(Rs2, Rp2), mf_set1(.5f)), 0.0f, 1.0f);
}

#if 0
static inline float fresnel(const float n1, const complex float n2, const float cosr)
{
  const complex float cost = csqrtf(1.0f - (n1/n2)*(n1/n2) * (1.0f - cosr*cosr));
  // fresnel for unpolarized light:
  const float Rs = cabsf((n1*cosr - n2*cost)/(n1*cosr + n2*cost));
  const float Rp = cabsf((n1*cost - n2*cosr)/(n1*cost + n2*cosr));
  return fminf(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}
#endif

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_reflect)) return mf_set1(0.0f);
  float wi[3];
  float wo[3];
  const float *n = p->v[v].hit.n;
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
  // these should always be positive:
  const float cos_in  = -dotproduct(n, wi);
  const float cos_out =  dotproduct(n, wo);
  if(cos_in < 0.0f) return mf_set1(0.0f);
  if(cos_out < 0.0f) return mf_set1(0.0f);
  float h[3];

  for(int k=0;k<3;k++) h[k] = wi[k] - wo[k];
  normalise(h);

  if(p->v[v].mode & s_specular)
  { // pure specular case reflect
    const float cosh = fabsf(dotproduct(h, n));
    if(cosh < HALFVEC_COS_THR) return mf_set1(0.0f);
    return mf_set1(1.0f);
  }

  float pdf = 1.0f / (4.0f * fabsf(dotproduct(wo, h)));

  pdf *= ggx_pdf_h(wi, h, n, p->v[v].shading.roughness);
  pdf /= fabsf(cos_out); // to projected solid angle measure
  // check for degenerate divisions etc:
  if(!(pdf > 0.0f)) return mf_set1(0.0f);
  return mf_set1(pdf);
}


mf_t sample(path_t *p, void* data)
{
  const metal_t *m = data;
  const int v = p->length-1; // current vertex.

  const float *n = p->v[v].hit.n;
  float ht[3] = {0.0, 0.0, 1.0}; // h in tangent space
  float pdf_h = 1.0f;
  float h[3] = {n[0], n[1], n[2]}; // init for smooth dielectric
  const float r = p->v[v].shading.roughness;
  if(r > GLOSSY_THR)
  {
    // flip incoming direction to point away from intersection point
    const float wit[3] = {
      - dotproduct(p->v[v].hit.a, p->e[v].omega),
      - dotproduct(p->v[v].hit.b, p->e[v].omega),
      - dotproduct(p->v[v].hit.n, p->e[v].omega)};
    ggx_sample_h(wit, r, r, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), ht);
    // world space:
    for(int k=0;k<3;k++) h[k] = ht[0]*p->v[v].hit.a[k] + ht[1]*p->v[v].hit.b[k] + ht[2]*n[k];
    pdf_h = ggx_pdf_h(p->e[v].omega, h, n, r);
  }
  float pdf = pdf_h;
  const float cosr = -dotproduct(p->e[v].omega, h);
  if(!(cosr > 0.0f)) return mf_set1(0.0f); // sampled a microfacet from the wrong side :(

  const mf_t n1 = p->e[v].vol.ior;
  mf_t n2, k2;
  fresnel_get_ior_mf(m->mat, p->lambda, &n2, &k2);
  const mf_t R = fresnel(n1, n2, k2, cosr);

  p->v[v].mode = s_reflect;
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k] + 2.0f * cosr * h[k];
  if(dotproduct(p->e[v+1].omega, n) <= 0.0f) return mf_set1(0.0f);
  pdf *= 1.0f/(4.0f * cosr);

  if(r > GLOSSY_THR)
  {
    p->v[v+1].pdf = mf_set1(pdf/fabsf(dotproduct(p->e[v+1].omega, n))); // convert to projected solid angle measure
    p->v[v].mode |= s_glossy;
    if(dotproduct(p->e[v+1].omega, n)*dotproduct(p->e[v+1].omega, h) < 0.0f) return mf_set1(0.0f);
    return mf_mul(R, mf_mul(p->v[v].shading.rg, mf_set1(ggx_shadowing_smith_G1(p->e[v+1].omega, n, h, p->v[v].shading.roughness))));
  }

  p->v[v].mode |= s_specular;
  return mf_mul(R, p->v[v].shading.rg);
}


mf_t brdf(path_t *p, int v, void *data)
{
  const metal_t *m = data;
  const float cos_in  = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  const float cos_out =  dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  const mf_t n1 = p->e[v].vol.ior;
  mf_t n2, k2;
  fresnel_get_ior_mf(m->mat, p->lambda, &n2, &k2);

  if(cos_out <= 0.0f || cos_in <= 0.0f) return mf_set1(0.0f);
  p->v[v].mode = s_reflect;

  const float r = p->v[v].shading.roughness;
  if(r > GLOSSY_THR)
    p->v[v].mode |= s_glossy;
  else
    p->v[v].mode |= s_specular;

  float h[3], cosh = 0.0f;
  for(int k=0;k<3;k++) h[k] = - p->e[v].omega[k] + p->e[v+1].omega[k];
  normalise(h);
  cosh = dotproduct(h, p->v[v].hit.n);
  if(cosh < 0.0f) return mf_set1(0.0f);

  // ggx half vector distribution:
  const float DG1 = ggx_pdf_h(p->e[v].omega, h, p->v[v].hit.n, p->v[v].shading.roughness);
  if(DG1 == 0) return mf_set1(0.0f);

  // fresnel on microfacet:
  const float cosr = -dotproduct(h, p->e[v].omega);
  if(cosr < 0.0f) return mf_set1(0.0f); // hit microfacet from back side
  // fresnel for unpolarized light:
  const mf_t R = fresnel(n1, n2, k2, cosr);
  const float G1 = ggx_shadowing_smith_G1(p->e[v+1].omega, p->v[v].hit.n, h, p->v[v].shading.roughness);

  if(cos_in == 0.0f) return mf_set1(0.0f);
  if(p->v[v].mode & s_glossy)
    return mf_mul(mf_mul(p->v[v].shading.rg, R), mf_set1(DG1 * G1 / (4.0f*fabsf(cosr*cos_out))));

  // now check dot to outgoing dir for pure specular case
  if(cosh < HALFVEC_COS_THR) return mf_set1(0.0f);
  return mf_mul(p->v[v].shading.rg, R);
}

