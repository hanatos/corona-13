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
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "colour.h"
#include "rgb2spec.h"
#include "mf.h"

// returns optional scale if colour was too bright
// depends on inited lut in global rt.rgb2spec
static inline float spectrum_rgb_to_coeff(const float rgb[3], float out[3])
{
  float col[3];
  float mul = MAX(MAX(rgb[0], rgb[1]), rgb[2]);
  if(mul == 0.0f || mul < 1.0f) mul = 1.0f;

  for(int k=0;k<3;k++) col[k] = rgb[k] / mul;
  rgb2spec_fetch(rt.rgb2spec, col, out);
  return mul;
}

static inline void spectrum_cauchy_from_abbe(const float n_d, const float V_d, float *A, float *B)
{   
  if(V_d == 0.0f)
  { 
    *A = n_d;
    *B = 0.0f;
    return;
  } 
  const float l_C = .6563f;
  const float l_F = .4861f;
  const float l_D = .587561f;
  const float c = (l_C*l_C * l_F*l_F)/(l_C*l_C - l_F*l_F);
  *B = (n_d - 1.0f)/V_d * c;
  *A = n_d - *B/(l_D*l_D);
} 

static inline mf_t spectrum_eta_from_abbe(const float n_d, const float V_d, const mf_t lambda)
{ 
  float A, B;
  spectrum_cauchy_from_abbe(n_d, V_d, &A, &B);

  // convert from micrometers to nanometers
  return mf_add(mf_set1(A), mf_div(mf_set1(B*1e6f), mf_mul(lambda, lambda)));
}

// 2-deg XYZ CMFs 1931
// http://cvrl.ioo.ucl.ac.uk/index.htm
static const int spectrum_sample_min = 360.0f;
static const int spectrum_sample_max = 830.0f;
static const int spectrum_xyz_lambda_min = 360;
static const int spectrum_xyz_lambda_max = 830;
static const int spectrum_xyz_step = 5;
// data + padded by one 0 so we don't need to clamp linear interpolation.
static const float spectrum_xyz_lut[] = {
0.000129900000,0.000003917000,0.000606100000,
0.000232100000,0.000006965000,0.001086000000,
0.000414900000,0.000012390000,0.001946000000,
0.000741600000,0.000022020000,0.003486000000,
0.001368000000,0.000039000000,0.006450001000,
0.002236000000,0.000064000000,0.010549990000,
0.004243000000,0.000120000000,0.020050010000,
0.007650000000,0.000217000000,0.036210000000,
0.014310000000,0.000396000000,0.067850010000,
0.023190000000,0.000640000000,0.110200000000,
0.043510000000,0.001210000000,0.207400000000,
0.077630000000,0.002180000000,0.371300000000,
0.134380000000,0.004000000000,0.645600000000,
0.214770000000,0.007300000000,1.039050100000,
0.283900000000,0.011600000000,1.385600000000,
0.328500000000,0.016840000000,1.622960000000,
0.348280000000,0.023000000000,1.747060000000,
0.348060000000,0.029800000000,1.782600000000,
0.336200000000,0.038000000000,1.772110000000,
0.318700000000,0.048000000000,1.744100000000,
0.290800000000,0.060000000000,1.669200000000,
0.251100000000,0.073900000000,1.528100000000,
0.195360000000,0.090980000000,1.287640000000,
0.142100000000,0.112600000000,1.041900000000,
0.095640000000,0.139020000000,0.812950100000,
0.057950010000,0.169300000000,0.616200000000,
0.032010000000,0.208020000000,0.465180000000,
0.014700000000,0.258600000000,0.353300000000,
0.004900000000,0.323000000000,0.272000000000,
0.002400000000,0.407300000000,0.212300000000,
0.009300000000,0.503000000000,0.158200000000,
0.029100000000,0.608200000000,0.111700000000,
0.063270000000,0.710000000000,0.078249990000,
0.109600000000,0.793200000000,0.057250010000,
0.165500000000,0.862000000000,0.042160000000,
0.225749900000,0.914850100000,0.029840000000,
0.290400000000,0.954000000000,0.020300000000,
0.359700000000,0.980300000000,0.013400000000,
0.433449900000,0.994950100000,0.008749999000,
0.512050100000,1.000000000000,0.005749999000,
0.594500000000,0.995000000000,0.003900000000,
0.678400000000,0.978600000000,0.002749999000,
0.762100000000,0.952000000000,0.002100000000,
0.842500000000,0.915400000000,0.001800000000,
0.916300000000,0.870000000000,0.001650001000,
0.978600000000,0.816300000000,0.001400000000,
1.026300000000,0.757000000000,0.001100000000,
1.056700000000,0.694900000000,0.001000000000,
1.062200000000,0.631000000000,0.000800000000,
1.045600000000,0.566800000000,0.000600000000,
1.002600000000,0.503000000000,0.000340000000,
0.938400000000,0.441200000000,0.000240000000,
0.854449900000,0.381000000000,0.000190000000,
0.751400000000,0.321000000000,0.000100000000,
0.642400000000,0.265000000000,0.000049999990,
0.541900000000,0.217000000000,0.000030000000,
0.447900000000,0.175000000000,0.000020000000,
0.360800000000,0.138200000000,0.000010000000,
0.283500000000,0.107000000000,0.000000000000,
0.218700000000,0.081600000000,0.000000000000,
0.164900000000,0.061000000000,0.000000000000,
0.121200000000,0.044580000000,0.000000000000,
0.087400000000,0.032000000000,0.000000000000,
0.063600000000,0.023200000000,0.000000000000,
0.046770000000,0.017000000000,0.000000000000,
0.032900000000,0.011920000000,0.000000000000,
0.022700000000,0.008210000000,0.000000000000,
0.015840000000,0.005723000000,0.000000000000,
0.011359160000,0.004102000000,0.000000000000,
0.008110916000,0.002929000000,0.000000000000,
0.005790346000,0.002091000000,0.000000000000,
0.004109457000,0.001484000000,0.000000000000,
0.002899327000,0.001047000000,0.000000000000,
0.002049190000,0.000740000000,0.000000000000,
0.001439971000,0.000520000000,0.000000000000,
0.000999949300,0.000361100000,0.000000000000,
0.000690078600,0.000249200000,0.000000000000,
0.000476021300,0.000171900000,0.000000000000,
0.000332301100,0.000120000000,0.000000000000,
0.000234826100,0.000084800000,0.000000000000,
0.000166150500,0.000060000000,0.000000000000,
0.000117413000,0.000042400000,0.000000000000,
0.000083075270,0.000030000000,0.000000000000,
0.000058706520,0.000021200000,0.000000000000,
0.000041509940,0.000014990000,0.000000000000,
0.000029353260,0.000010600000,0.000000000000,
0.000020673830,0.000007465700,0.000000000000,
0.000014559770,0.000005257800,0.000000000000,
0.000010253980,0.000003702900,0.000000000000,
0.000007221456,0.000002607800,0.000000000000,
0.000005085868,0.000001836600,0.000000000000,
0.000003581652,0.000001293400,0.000000000000,
0.000002522525,0.000000910930,0.000000000000,
0.000001776509,0.000000641530,0.000000000000,
0.000001251141,0.000000451810,0.000000000000,
0,0,0,
};

static inline void spectrum_xyz(const float lambda, float *xyz)
{
  assert(lambda >= spectrum_xyz_lambda_min);
  assert(lambda <= spectrum_xyz_lambda_max); // we allow == here since the cmf are padded with one additional 0 entry to catch cases where the random number == 1.0.
  //if(lambda < spectrum_xyz_lambda_min || lambda >= spectrum_xyz_lambda_max) { xyz[0] = xyz[1] = xyz[2] = 0.0f; return; }
  float f = (lambda - spectrum_xyz_lambda_min)/spectrum_xyz_step;
  int i = (int)f;
  f -= i;
  xyz[0] = (1-f)*spectrum_xyz_lut[3*i+0] + f*spectrum_xyz_lut[3*(i+1)+0];
  xyz[1] = (1-f)*spectrum_xyz_lut[3*i+1] + f*spectrum_xyz_lut[3*(i+1)+1];
  xyz[2] = (1-f)*spectrum_xyz_lut[3*i+2] + f*spectrum_xyz_lut[3*(i+1)+2];
}

static inline void spectrum_p_to_xyz(const mf_t lambda, const mf_t p, float *xyz)
{
  float b[3];
  const float *pf = (const float *)&p;
  const float *lf = (const float *)&lambda;
  xyz[0] = xyz[1] = xyz[2] = 0.0f;
  for(int l=0;l<MF_COUNT;l++)
  {
    spectrum_xyz(lf[l], b);
    for(int k=0;k<3;k++) xyz[k] += b[k]*pf[l];
  }
}

static inline void spectrum_p_to_camera(const mf_t lambda, const mf_t p, float *cam)
{
  float xyz[3];
  spectrum_p_to_xyz(lambda, p, xyz);
  colour_xyz_to_camera(xyz, cam);
}

// p = 1/(720-380)
static inline mf_t spectrum_sample_lambda(const mf_t rand, mf_t *pdf)
{
  if(pdf) *pdf = mf_set1(1.0f/(spectrum_sample_max - spectrum_sample_min));
  return mf_add(mf_set1(spectrum_sample_min), mf_mul(mf_set1(spectrum_sample_max - spectrum_sample_min), rand));
}

static inline mf_t spectrum_lambda_pdf(const mf_t lambda)
{
  return mf_set1(1.0f/(spectrum_sample_max - spectrum_sample_min));
}

// #define SPECTRUM_LAMBDA_STEP 5.0f
#define SPECTRUM_LAMBDA_STEP 50.0f
static inline float spectrum_mutate(const float lambda, float rand, float *pdf)
{
  const float step = SPECTRUM_LAMBDA_STEP; // in nanometers
  if(rand > 0.5f) rand = - 2.0f*step*(rand - .5f);
  else rand = rand * 2.0f * step;
  // mirror on boundary
  float l = lambda + rand;
  if(l < spectrum_sample_min) l = 2.0f*spectrum_sample_min - l;
  if(l > spectrum_sample_max) l = 2.0f*spectrum_sample_max - l;
  if(pdf) *pdf = .5f/step;
  return l;
}

static inline float spectrum_pdf_mutate(const float l1, const float l2)
{
  const float step = SPECTRUM_LAMBDA_STEP;
  const float e = 1e-3f;
  if(fabsf(l1 - l2) <= step+e ||
     fabsf(l1-spectrum_sample_max) + fabsf(l2-spectrum_sample_min) <= step + e ||
     fabsf(l2-spectrum_sample_max) + fabsf(l1-spectrum_sample_min) <= step + e)
    return .5f/step;
  return 0.0f;
}
#undef SPECTRUM_LAMBDA_STEP
