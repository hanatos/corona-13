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
#ifndef SAMPLER_COMMON_H
#define SAMPLER_COMMON_H

#include "corona_common.h"
#include <assert.h>

// some useful Monte Carlo and quasi-Monte Carlo functions

static inline float sample_mutate_rand(float x, const float rand, const float amount)
{
  // FIXME: this is not symmetric by one ulp
  const float dx = (2.f*rand - 1.f)*amount;
  float x1 = x + dx;
  if(dx < 0.0f)
    return (x1 < 0) ? x1 + 1 : x1;
  else
    return (x1 > 1) ? x1 - 1 : x1;
}

#ifndef __cplusplus // :(
#include <complex.h>
// transforms a 2d point on the disk such that (0,0) maps to the given pivot.
// this is an automorphism (i.e. never leaves the disk).
static inline void moebius_sample(
    const float r1,  // input random numbers
    const float r2,
    float x0,        // pivot of the transform
    float y0,
    float k_u,       // exponent
    float k_v,       // exponent
    float *u,        // output point
    float *v)
{
  // const float cos_phi = cosf(r1*M_PI*2.0f), sin_phi = sinf(r1*M_PI*2.0f);
  // const float cos_theta = powf(1.0f - r2, 1.0f/(k+1.0f));
  // const float sin_theta = sqrtf(fmaxf(0.0f, 1.0f - cos_theta*cos_theta));
  float cos_phi = sqrtf(k_v+1.0f)*cosf(r1*M_PI*2.0f), sin_phi = sqrtf(k_u+1.0f)*sinf(r1*M_PI*2.0f);
  const float ilen = 1.0f/sqrtf(cos_phi*cos_phi + sin_phi*sin_phi);
  cos_phi *= ilen; sin_phi *= ilen;
  const float cos_theta = powf(1.0f - r2, 1.0f/(k_u * cos_phi*cos_phi + k_v * sin_phi*sin_phi));
  const float sin_theta = sqrtf(fmaxf(0.0f, 1.0f - cos_theta*cos_theta));

  const float complex z = cos_phi*sin_theta + I*sin_phi*sin_theta;
  const float complex z0 = x0 + I*y0;

  const float complex den = (conjf(z0)*z - 1);
  if(fabsf(crealf(den)) < 1e-8 && fabsf(cimagf(den)) < 1e-8)
  { // arrgg
    *u = *v = 0.0f;
    return;
  }
  // moebius transform w = f(z):
  const float complex w = (z - z0)/den;

  // store in outgoing direction:
  *u = crealf(w);
  *v = cimagf(w);
}

static inline float moebius_eval(
    float u,   // previously transformed point
    float v,
    float x0,  // pivot point
    float y0,
    float k_u, // exponent
    float k_v)
{
  const float complex z0 = x0 + I*y0; // pivot
  const float complex w  = u + I*v;   // input point
  const float complex z = (z0 - w)/(1.0f - conjf(z0)*w); // backtransform before moebius
  const float complex den = (conjf(z0)*z0 - 1);
  if(fabsf(crealf(den)) < 1e-8 && fabsf(cimagf(den)) < 1e-8)
    return 0.0f;
  const float complex dfdx = (z*conjf(z0) - 1.0f)*(z*conjf(z0) - 1.0f)/den;
  const float cos_theta = sqrtf(fmaxf(0.0f, 1.0f - crealf(z*conjf(z))));
  // return powf(cos_theta, k) * (k+1.0f)/(2.0f*M_PI) * crealf(conjf(dfdx)*dfdx);
  const float len2 = crealf(conjf(z)*z);
  float cos2_phi = 1.0f, sin2_phi = 0.0f;
  if(len2 > 1e-8f)
  {
    cos2_phi = crealf(z)*crealf(z)/len2;
    sin2_phi = cimagf(z)*cimagf(z)/len2;
  }
  return powf(cos_theta, cos2_phi * k_u + sin2_phi * k_v) * sqrtf((k_u+1.0f)*(k_v+1.0f))/(2.0f*M_PI) * crealf(conjf(dfdx)*dfdx);
}
#endif

static inline float sample_cubic_bspline(
    const float r1, const float r2,
    const float r3, const float r4)
{
  return -4.0f/2.0f + r1 + r2 + r3 + r4;
}

static inline float sample_cubic_bspline_pdf(
    const float x)
{
  const float ax = fabsf(x);
  if(ax >= 2) return 0.0f;
  if(ax >= 1) return (2.0f-ax)*(2.0f-ax)*(2.0f-ax)/6.0f;
  return (4.0f - 6.0f*ax*ax + 3.0f*ax*ax*ax)/6.0f;
}

// samples two gaussian distributed numbers
static inline void sample_gaussian(
    const float r1,
    const float r2,
    float *g1,
    float *g2)
{
  const float r = sqrtf(fmaxf(0.0f, -2.0f*logf(1.0f-r1)));
  float sint, cost;
  common_sincosf(2.0f*M_PI*r2, &sint, &cost);
  *g1 = sint * r;
  *g2 = cost * r;
}


// uniformly sample the unit sphere, p = 1/4pi
static inline void sample_sphere(float *x, float *y, float *z, const float x1, const float x2)
{
  *z = 1.f - 2.f*x1;
  const float r = sqrtf(1.f - *z**z);
  const float phi = 2.f*M_PI*x2;
  *x = r * cosf(phi);
  *y = r * sinf(phi);
}

// sample hemisphere uniformly, p = 1/2pi
static inline void sample_hemisphere(float *x, float *y, float *z, const float x1, const float x2)
{
  *z = 1.f - x1;
  const float r = sqrtf(1.f - *z**z);
  const float phi = 2.f*M_PI*x2;
  *x = r * cosf(phi);
  *y = r * sinf(phi);
}

// sample hemisphere, cos lobe, p = cos(theta)/pi
static inline void sample_cos(float *x, float *y, float *z, const float x1, const float x2)
{
  const float su = sqrtf(x1);
  *x = su*cosf(2.f*M_PI*x2);
  *y = su*sinf(2.f*M_PI*x2);
  *z = sqrtf(1.0 - x1);
}

// sample hemisphere, cos^k lobe, p = cos^k(theta) (k+1)/2pi
static inline void sample_cos_k(float *x, float *y, float *z, const float k, const float x1, const float x2)
{
  const float r1 = x1 * 2.0f * M_PI;
  const float cos_theta = powf(1.0f - x2, 1.0f/(k+1));
  const float sin_theta = sqrtf(MAX(0.0f, 1.0f - cos_theta*cos_theta));
  *x = cosf(r1) * sin_theta;
  *y = sinf(r1) * sin_theta;
  *z = cos_theta;
}

// fibonacci rank-1 lattice
static inline void sample_rank1_fib21(const unsigned int i, float *x, float *y)
{
  *x = i/21.0f;
  *y = i*(13.0f/21.0f);
  *x -= (int)*x;
  *y -= (int)*y;
}

static inline uint32_t sample_cdfd(const double *cdf, const int num, const float rand)
{
  // unsigned int t = 0;
  // for(;cdf[t] < rand && t < num;t++);
  unsigned int min = 0, max = num;
  unsigned int t = max/2;
  while (t != min)
  {
    if(cdf[t] <= rand) min = t;
    else max = t;
    t = (min + max)/2;
  }
  // last step: decide between min and max one more time (min is rounding default),
  // but consider that if max is still out of bounds, it's invalid.
  // (rand == 1.0 and a cdf[0]=1, num=1 would break otherwise)
  if(max < num && cdf[t] <= rand) t = max;
  assert(t==0 || cdf[t-1] <= rand);
  assert(max == num || cdf[t] > rand); // max == num special case is for above mentioned edge case, too
  assert(t<num);
  return t;
}

static inline uint32_t sample_cdf(const float *cdf, const int num, const float rand)
{
  // unsigned int t = 0;
  // for(;cdf[t] < rand && t < num;t++);
  unsigned int min = 0, max = num;
  unsigned int t = max/2;
  while (t != min)
  {
    if(cdf[t] <= rand) min = t;
    else max = t;
    t = (min + max)/2;
  }
  // last step: decide between min and max one more time (min is rounding default),
  // but consider that if max is still out of bounds, it's invalid.
  // (rand == 1.0 and a cdf[0]=1, num=1 would break otherwise)
  if(max < num && cdf[t] <= rand) t = max;
  assert(t==0 || cdf[t-1] <= rand);
  assert(max == num || cdf[t] > rand); // max == num special case is for above mentioned edge case, too
  assert(t<num);
  return t;
}

// sample a henyey greenstein lobe with mean cosine g
static inline void sample_hg(
    const float g,  // mean cosine
    const float r1, // one uniform [0,1) random number
    const float r2, // another of those
    float *out,     // output in local tangent space coords
    float *pdf)     // if != 0, computes and returns pdf here.
{
  if(g == 0.0f)
  {
    sample_sphere(out, out+1, out+2, r1, r2);
    if(pdf) *pdf = 1.0f/(4.0f*M_PI);
    return;
  }
  const float sqr = (1.0f-g*g)/(1.0f-g*(2.0f*r1-1));
  const float cos_theta = 1.0f/(2.0f*g) * (1.0f + g*g - sqr*sqr);
  const float phi = 2.0f*M_PI*r2;
  const float l = sqrtf(fmaxf(0.0f, 1.0f-cos_theta*cos_theta));
  out[0] = cos_theta;
  out[1] = cosf(phi)*l;
  out[2] = sinf(phi)*l;
  normalise(out);
  if(pdf) *pdf = 1.0f/(4.0f*M_PI) * (1.0f-g)/powf(1.0f + g*g - 2.0f*g*cos_theta, 3.0f/2.0f);
}

static inline void sample_inverse_hg(
    const float g,
    const float *wo,
    float *r1,
    float *r2)
{
  if(g == 0.0f)
  {
    *r1 = (1.0f-wo[2])/2.0f;
    const float phi = atan2f(wo[1], wo[0]);
    *r2 = phi / (2.0f*M_PI);
    return;
  }
  const float cos_theta = wo[0];
  const float phi = atan2f(wo[2], wo[1]);
  *r2 = phi / (2.0f*M_PI);
  *r1 = (1.0f-g*g)/(2.0f*g)*(1.0f/sqrtf(1.0f+g*g-2.0f*g*cos_theta) - 1.0f/(1.0f+g));
}

static inline float sample_eval_hg(
    const float g,
    const float *wi,
    const float *wo)
{
  if(g == 0.0f) return 1.0f/(4.0f*M_PI);
  const float cos_theta = dotproduct(wi, wo);
  return 1.0f/(4.0f*M_PI) * (1.0f-g)/powf(1.0f + g*g - 2.0f*g*cos_theta, 3.0f/2.0f);
}



#endif
