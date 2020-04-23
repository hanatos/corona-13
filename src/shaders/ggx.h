#pragma once
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

#include "mf.h"

// eric and eugene's `importance sampling microfacet based bsdfs using the distribution of visible normals',
// ggx + smith shadowing variant.

// TODO: copy some of the anisotropic bits in the microfacet header for
// TODO: multiple scattering here, too

// unidirectional shadowing, input in world space
static inline float ggx_shadowing_smith_G1(const float *w, const float *n, const float *h, const float roughness)
{
  const float r2 = roughness * roughness;
  const float cos_th = fabsf(dotproduct(w, n));
  const float sin_th = sqrtf(fmaxf(0.0f, 1.0f-cos_th*cos_th));
  const float tan_th = sin_th/cos_th;
  return 2.0f/(1.0f + sqrtf(1.0f + r2 * tan_th * tan_th));
}

static inline mf_t ggx_shadowing_smith_G1_mf(
    const mf_t cos_wn,
    const float roughness)
{
  const mf_t r2 = mf_set1(roughness * roughness);
  const mf_t sin_wn = mf_sqrt(mf_clamp(mf_sub(mf_set1(1.0f), mf_mul(cos_wn, cos_wn)), 0.0f, 1.0f));
  const mf_t tan_th = mf_div(sin_wn, cos_wn);
  return mf_div(mf_set1(2.0f), mf_add(mf_set1(1.0f), mf_sqrt(mf_add(mf_set1(1.0f), mf_mul(mf_mul(r2,  tan_th), tan_th)))));
}

// bidirectional shadowing, input in world space
static inline float ggx_shadowing_smith_G2(const float *wi, const float *wo, const float *n, const float *h, const float roughness)
{
  // check sidedness of micro vs macro normal:
  if(dotproduct(wi, n)*dotproduct(wi, h) < 0.0f) return 0.0f;
  if(dotproduct(wo, n)*dotproduct(wo, h) < 0.0f) return 0.0f;
  return ggx_shadowing_smith_G1(wi, n, h, roughness) * ggx_shadowing_smith_G1(wo, n, h, roughness);
}


// helper to stretch pdf
static inline void _ggx_sample11(
    // input
    const float tan_theta_i, // tangent of theta incident tan(normal, wi)
    float U1, float U2, // random numbers
    // output
    float *slope_x,
    float *slope_y)
{
  // special case (normal incidence)
  // XXX FIXME: this threshold is unwisely chosen, it causes inconsistency between sampling and bsdf for normal incidence!
  if(tan_theta_i < 0.0001f)
  {
    const float r = sqrtf(U1/fmaxf(1e-8f, 1-U1));
    const float phi = 2.0f * M_PI * U2;
    *slope_x = r * cosf(phi);
    *slope_y = r * sinf(phi);
    return;
  }
  // precomputations
  const float a = 1.0f / tan_theta_i;
  const float G1 = 2.0f / (1.0f + sqrtf(1.0f+1.0f/(a*a)));
  // sample slope_x
  const float A = 2.0f*U1/G1 - 1.0f;
  const float tmp = 1.0f / (A*A-1.0f);
  const float B = tan_theta_i;
  const float D = sqrtf(fmaxf(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
  float slope_x_1 = B*tmp - D;
  float slope_x_2 = B*tmp + D;
  // fix cases where B = inf and tmp = 0 and tmp = inf etc.
  if(!(fabsf(slope_x_1) < FLT_MAX)) slope_x_1 = 0.0f;
  if(!(fabsf(slope_x_2) < FLT_MAX)) slope_x_2 = 0.0f;
  *slope_x = (A < 0.0f || slope_x_2*tan_theta_i > 1.0f) ? slope_x_1 : slope_x_2;
  // sample slope_y
  float S;
  if(U2 > 0.5f)
  {
    S = 1.0f;
    U2 = 2.0f*(U2-0.5f);
  }
  else
  {
    S = -1.0f;
    U2 = 2.0f*(0.5f-U2);
  }
  // improved fit from mitsuba:
  const float z = (U2 * (U2 * (U2 * (-0.365728915865723f) + 0.790235037209296f) -
        0.424965825137544f) + 0.000152998850436920f) /
    (U2 * (U2 * (U2 * (U2 * 0.169507819808272f - 0.397203533833404f) -
                 0.232500544458471f) + 1.0f) - 0.539825872510702f);
  *slope_y = S * z * sqrtf(1.0+*slope_x * *slope_x);
  assert(*slope_y == *slope_y);
}

// sample ggx microfacet normal
// remaining Monte Carlo weight is:
// ggx_shadowing_smith_G1(wo, n, h, roughness);
static inline void ggx_sample_h(
    // input
    const float wi[3], // incident direction, tangent space of hit point, pointing away from surface
    const float roughness_x, const float roughness_y, // anisotropic roughness
    const float U1, const float U2, // random numbers
    // output
    float h[3]) // micronormal in tangent space of hit point
{
  // 1. stretch omega_i
  float wi_[3];
  wi_[0] = roughness_x * wi[0];
  wi_[1] = roughness_y * wi[1];
  wi_[2] = fabsf(wi[2]);
  normalise(wi_);
  // get polar coordinates of omega_i_
  float tan_theta = 0.0f;
  float sin_phi = 0.0f, cos_phi = 1.0f;
  if (wi_[2] < 0.99999)
  {
    const float len = sqrtf(wi_[0]*wi_[0] + wi_[1]*wi_[1]);
    tan_theta = len/wi_[2];
    sin_phi = wi_[1]/len;
    cos_phi = wi_[0]/len;
  }
  // 2. sample P22_{wi}(x_slope, y_slope, 1, 1)
  float slope_x, slope_y;
  _ggx_sample11(
      tan_theta,
      U1, U2,
      &slope_x, &slope_y);
  // 3. rotate
  float tmp = cos_phi*slope_x - sin_phi*slope_y;
  slope_y = sin_phi*slope_x + cos_phi*slope_y;
  slope_x = tmp;
  // 4. unstretch
  slope_x = roughness_x * slope_x;
  slope_y = roughness_y * slope_y;
  // 5. compute normal
  float inv_h = sqrtf(slope_x*slope_x + slope_y*slope_y + 1.0);
  h[0] = -slope_x/inv_h;
  h[1] = -slope_y/inv_h;
  h[2] = 1.0/inv_h;
  if(!(inv_h > 0.0))
  {
    h[0] = h[2] = 0.0f;
    h[1] = 1.0f;
  }
}

// visible micro normal distribution in solid angle measure.
// eq (2) in Heitz/d'Eon paper:
// D_{\omega_i}(\omega_m) =
static inline float ggx_pdf_h(
    const float *wi,  // omega in, world space
    const float *h,   // half vector (micro facet normal), world space
    const float *n,   // normal, world space
    const float roughness)
{
  const float r2 = roughness*roughness;
  const float cos_th = fabsf(dotproduct(h, n));
  const float sin_th = sqrtf(fmaxf(0.0f, 1.0f-cos_th*cos_th));
  const float tan_th = sin_th/cos_th; // tan theta h
  const float D_h = r2/(M_PI * cos_th*cos_th*cos_th*cos_th * (r2 + tan_th*tan_th)*(r2 + tan_th*tan_th));
  const float G1 = ggx_shadowing_smith_G1(wi, n, h, roughness);

  // distribution of visible normals:
  return fabsf(G1 * dotproduct(wi, h) * D_h / dotproduct(wi, n));
}

static inline mf_t ggx_pdf_h_mf(
    const mf_t cosh,   // = dot(n, h)
    const mf_t cos_in, // = dot(wi, n)
    const mf_t cosr,   // = dot(wi, h)
    const float roughness)
{
  const mf_t r2 = mf_set1(roughness*roughness);
  const mf_t cosh2 = mf_mul(cosh, cosh);
  const mf_t sin_th = mf_sqrt(mf_clamp(mf_sub(mf_set1(1.0f), cosh2), 0.0f, 1.0f));
  const mf_t tan_th = mf_div(sin_th, mf_abs(cosh)); // tan theta h
  const mf_t den = mf_fma(tan_th, tan_th, r2);
  const mf_t ct4 = mf_mul(cosh2, cosh2);
  const mf_t D_h = mf_div(r2, mf_mul(mf_mul(mf_set1(M_PI), ct4), mf_mul(den, den)));
  const mf_t G1 = ggx_shadowing_smith_G1_mf(cos_in, roughness);

  // distribution of visible normals:
  return mf_abs(mf_mul(mf_mul(G1, cosr), mf_div(D_h, cos_in)));
}

static inline int ggx_inverse_sample(
    // input
    const float wi[3], // incident direction, tangent space of hit point, pointing away from surface
    const float roughness_x, const float roughness_y, // anisotropic roughness
    const float h[3], // micronormal in tangent space of hit point
    // output
    float *U1, float *U2) // random numbers
{
  // 1. stretch omega_i
  float wi_[3];
  wi_[0] = roughness_x * wi[0];
  wi_[1] = roughness_y * wi[1];
  wi_[2] = fabsf(wi[2]);
  normalise(wi_);
  // get polar coordinates of omega_i_
  float tan_theta_i = 0.0f;
  float sin_phi = 0.0f, cos_phi = 1.0f;
  float sin_theta_i = sqrtf(wi_[0]*wi_[0] + wi_[1]*wi_[1]);
  float cos_theta_i = wi_[2];
  if (wi_[2] < 0.99999)
  {
    tan_theta_i = sin_theta_i/cos_theta_i;
    sin_phi = wi_[1]/sin_theta_i;
    cos_phi = wi_[0]/sin_theta_i;
  }
  // 5. compute slopes from normal
  assert(h[2] > 0.0f); // in fact,
  assert(dotproduct(wi, h) > 0.0f);
  float slope_x = -h[0]/h[2];
  float slope_y = -h[1]/h[2];
  // 4. stretch
  slope_x = slope_x / roughness_x;
  slope_y = slope_y / roughness_y;
  // 3. rotate backwards
  float tmp = cos_phi*slope_x + sin_phi*slope_y;
  slope_y = -sin_phi*slope_x + cos_phi*slope_y;
  slope_x = tmp;

  // 2. invert the sampling of P22_{wi}(x_slope, y_slope, 1, 1)
  const float z = slope_y / sqrtf(1.0+slope_x*slope_x);
  *U2 = (((z*z+1.0) * atanf(z) + z)/(2.0*z*z+2.0) + M_PI/4.0)/(M_PI/2.0);
  const float G1 = 2.0f / (1.0f + sqrtf(1.0f+tan_theta_i*tan_theta_i));
  *U1 = G1/(2.0*cos_theta_i) *
    ((sin_theta_i + slope_x * cos_theta_i)/sqrtf(1.0+slope_x*slope_x) + cos_theta_i);
  return 0;
}

