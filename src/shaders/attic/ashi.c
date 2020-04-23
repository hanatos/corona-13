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
#include <assert.h>


float powf5(const float x)
{
  const float x2 = x*x;
  return x2*x*x2;
}

static inline float get_fresnel(const float rs, const float cos_theta_in, const float cos_theta_out)
{
  return fminf(1.0f, rs + rs*powf5((1.0f - .5f*cos_theta_in)*(1.0f - .5f*cos_theta_out)));
}

float prepare(path_t *p, int v, void *data)
{
  p->v[v].material_modes = s_reflect | s_glossy;
  return 1.0f;
}

static inline float ashi_pdf(path_t *p, int e1, int v, int e2, void *data)
{
  const float nu = (2.0f/p->v[v].shading.roughness - 2.0f);
  float wi[3];
  float wo[3];
  float n[3];
  // TODO: factor out to shader.c. we all have to do this :(
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
      n[k] = p->v[v].hit.n[k];
    }
  }

  float doto = dotproduct(n, wo);
  float doti = - dotproduct(n, wi);
  if(doti * doto <= 0) return 0.0f;
  doto = fabsf(doto);
  doti = fabsf(doti);

  const float p_s = p->v[v].shading.rs + (1.0f-p->v[v].shading.rs)*powf5(1.0f-doti);
  float h[3];
  for(int k=0;k<3;k++) h[k] = wo[k] - wi[k];
  normalise(h);
  //const float *u = hit->a, *v = hit->b;
  // const float cos = dotproduct(h, u);
  // const float sin = dotproduct(h, v);
  // return p_s*sqrtf((nu+1)*(nv+1))*powf(dotproduct(n, h), nu*cos*cos + nv*sin*sin)*fmaxf(doto, doti)/(-8.0f*M_PI*dotproduct(omega_in, h)) + (1.0f-p_s)/M_PI;
  // isotropic case:
  // formula as in ashikhmin's paper /doto to convert to projected solid angle
  return p_s*(nu+1)*powf(fabsf(dotproduct(n, h)), nu)/(8.0f*M_PI*doto*fabsf(dotproduct(wi, h))) + (1.0-p_s)/M_PI;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  return ashi_pdf(p, e1, v, e2, data);
}

float sample(path_t *p, void *data)
{
  const int vi = p->length-1;
  const float nu = (2.0f/p->v[vi].shading.roughness - 2.0f);
  // float nv = *((float*)data + 1);
  float *n = p->v[vi].hit.n;
  float *u = p->v[vi].hit.a, *v = p->v[vi].hit.b;
  float *wi = p->e[p->length-1].omega;
  float *wo = p->e[p->length].omega;
  const float rs = p->v[vi].shading.rs;
  const float rd = p->v[vi].shading.rd;

  const float x1 = pointsampler(p, s_dim_omega_x);
  const float x2 = pointsampler(p, s_dim_omega_y);

  float nwi = - n[0]*wi[0] - n[1]*wi[1] - n[2]*wi[2];
  // importance sample by fresnel term (not reciprocal!)
  const float p_s = rs + (1.0f-rs)*powf5(1.0f-nwi);
  //const float p_d = (1.0f - p_s)* hit->rd * (1.0f-hit->rs);
  float throughput = 0.0f;
  if(pointsampler(p, s_dim_scatter_mode) < p_s)
  {
    // isotropic specular part
    float phi = 2.0f*M_PI*x1;
#if 1
    float cos_theta = powf(1 - x2, 1.0f/(nu + 1.0f));
    float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
    const float cos_phi = cosf(phi), sin_phi = sinf(phi);
#else // aniso:
    float cos_phi = sqrtf(nv+1.0f)*cosf(phi), sin_phi = sqrtf(nu+1.0f)*sinf(phi);
    const float ilen = 1.0f/sqrtf(cos_phi*cos_phi + sin_phi*sin_phi);
    cos_phi *= ilen; sin_phi *= ilen;
#endif

    float h[3];
    for(int k=0; k<3; k++) h[k] = sin_theta*cos_phi*u[k] + sin_theta*sin_phi*v[k] + cos_theta*n[k];

    float hw = - h[0]*wi[0] - h[1]*wi[1] - h[2]*wi[2];
    for(int k=0;k<3;k++) wo[k] = wi[k] + 2*hw*h[k];
    const float nwo = n[0]*wo[0] + n[1]*wo[1] + n[2]*wo[2];

    if(!(nwo > 0.0f)) return 0.0f;
    // p->v[vi].pdf *= p_s*(nu+1)*powf(cos_theta, nu)/(8.0f*M_PI*nwo*fabsf(hw));
    throughput = get_fresnel(rs, nwi, nwo) * nwo / (p_s * fmaxf(nwi, nwo));
  }
  else //if(rr < p_s + p_d)
  {
    // cosine distribution:
    float phi = 2 * M_PI * x1;
    for(int k=0; k<3; k++) wo[k] = sqrtf(x2)*cosf(phi)*u[k] + sqrtf(x2)*sinf(phi)*v[k] + sqrtf(1.0f-x2)*n[k];
    const float nwo = n[0]*wo[0] + n[1]*wo[1] + n[2]*wo[2];
    if(!(nwo > 0.0f)) return 0.0f;
#if 0
    float h[3];
    for(int k=0;k<3;k++) h[k] = wo[k] - wi[k];
    normalise(h);
    const float cos_theta = dotproduct(h, n);
    float hw =  - h[0]*wi[0] - h[1]*wi[1] - h[2]*wi[2];
    p->v[vi].pdf *= p_s*(nu+1)*powf(cos_theta, nu)/(8.0f*M_PI*nwo*fabsf(hw)) + (1.0-p_s)/M_PI;
#endif
    throughput = ((1.0f - get_fresnel(rs, nwi, nwo))/(1.0f- p_s)) * rd * (1.0f-rs)*(28.0f/23.0f);
  }

  p->v[vi].mode |= s_glossy | s_reflect;
  p->v[vi+1].pdf = ashi_pdf(p, vi, vi, vi+1, data);
  return throughput;
}

float brdf(path_t *p, int v, void *data)
{
  const float nu = (2.0f/p->v[v].shading.roughness - 2.0f);
  // const float nv = *((float*)data + 1);
  float *n = p->v[v].hit.n;
  const float rs = p->v[v].shading.rs;
  const float rd = p->v[v].shading.rd;
  float *k1 = p->e[v].omega;
  float *k2 = p->e[v+1].omega;
  float nk1 = n[0]*k1[0] + n[1]*k1[1] + n[2]*k1[2];
  float nk2 = n[0]*k2[0] + n[1]*k2[1] + n[2]*k2[2];
  if(nk1*nk2 >= 0.0f) return 0.0f;
  nk1 = fabsf(nk1);
  nk2 = fabsf(nk2);
  float h[3] = {k2[0]-k1[0], k2[1]-k1[1], k2[2]-k1[2]};
  float hlen = sqrtf(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
  for(int k=0;k<3;k++) h[k] /= hlen;
  float hk = fabsf(h[0]*k2[0] + h[1]*k2[1] + h[2]*k2[2]);
  const float fresnel = get_fresnel(rs, nk1, nk2);
  //diffuse
  const float pdconst = (28.0f/(23.0f*M_PI)) * (1.0f - fresnel);
  float brdf = rd * (1.0f-rs) * pdconst;
  p->v[v].mode |= s_reflect | s_glossy;

#if 1
  float hn = fmaxf(0.0f, h[0]*n[0] + h[1]*n[1] + h[2]*n[2]);
  const float psconst = (nu+1)*powf(hn, nu) / (8*M_PI*fmaxf(nk1,nk2)*hk);
#else // aniso
  float *u = p->v[v].hit.a, *v = p->v[v].hit.b;
  float hu = fabsf(h[0]*u[0] + h[1]*u[1] + h[2]*u[2]);
  float hv = fabsf(h[0]*v[0] + h[1]*v[1] + h[2]*v[2]);
  const float psconst = sqrtf((nu+1)*(nv+1))*powf(hn,(nu*hu*hu + nv*hv*hv)/(1-hn*hn))
           / (8*M_PI*fmaxf(nk1,nk2)*hk);
#endif
  brdf += psconst * fresnel;
  return brdf;
}

extern int init(FILE *f, void **data)
{
  *data = 0;
  int dreggn = fscanf(f, "%*[^\n]\n");
  return dreggn == -1;
}


