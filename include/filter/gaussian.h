/*
    This file is part of corona-13.

    copyright (c) 2018 johannes hanika.

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

#include "corona_common.h"
#include "matrix2.h"
#include "fakegaussian.h"
#include "svd2.h"

static inline void filter_gaussian_splat(
    framebuffer_t *fb,
    const float i,
    const float j,
    const float *mu,
    const float *S,
    const float *col)
{
  const int wd = fb->header->width, ht = fb->header->height;
#if 1
  // use svd and regularise eigenvalues:
  double ev_d[2] = {0};
  double Et_d[4] = {0};
  double S_d[4] = { S[0], S[2], S[1], S[3] };
  getEVDSymmetric2x2(ev_d, Et_d, S_d);
  float ev[] = { ev_d[0], ev_d[1] };
  float Et[] = {Et_d[0], Et_d[1], Et_d[2], Et_d[3]};
  float B[4] = {0};
  float D[4] = {0};
  // clamp eigenvalues
  ev[0] = CLAMP(sqrtf(ev[0]), 0.5f, 7.0f);
  ev[1] = CLAMP(sqrtf(ev[1]), 0.5f, 7.0f);
  D[0] = 1.f/ev[0];
  D[3] = 1.f/ev[1];
  mat2_mul(D, Et, B);
  const float det = 1.0f/fabsf(ev[0]*ev[1]);
#else
  // cholesky decomposition of covariance matrix:
  float Si[4];
  mat2_invert(S, Si);
  float B[4];
  B[0] = sqrtf(MAX(0.001f, Si[0]));
  B[1] = Si[1] / B[0];
  B[2] = 0.0f;
  B[3] = sqrtf(MAX(0.001f, Si[3] - B[1]*B[1]));

  B[0] = MAX(0.5, B[0]);
  B[3] = MAX(0.5, B[3]);
  // actually need decomp of inverse cov matrix:
  // mat2_invert(B, Bi);
  float det = fabsf(B[0]*B[3]);
  float ev[2] = {1.0f, 1.0f};
  float *Et = B;
#endif

  // bound the anisotropic gaussian:
  float bound[4] = {0, 0, 0, 0};
  float b0[] = {-4, 4}, b1[] = {4, 4};
  float tb0[2], tb1[2];
  for(int k=0;k<2;k++)
  {
    b0[k] *= ev[k];
    b1[k] *= ev[k];
  }
  float E[4];
  mat2_transpose(Et, E);
  mat2_mulv(E, b0, tb0);
  mat2_mulv(E, b1, tb1);
  for(int k=0;k<2;k++)
  {
    if(bound[k]   > tb0[k]) bound[k]   = tb0[k];
    if(bound[3+k] < tb0[k]) bound[3+k] = tb0[k];
    if(bound[k]   > tb1[k]) bound[k]   = tb1[k];
    if(bound[3+k] < tb1[k]) bound[3+k] = tb1[k];
  }

  // apply symmetry:
  bound[0] = MIN(bound[0], -bound[2]);
  bound[1] = MIN(bound[1], -bound[3]);
  bound[2] = MAX(bound[2], -bound[0]);
  bound[3] = MAX(bound[3], -bound[1]);
  // some last resort clamping
  float w = MIN(20, (bound[2]-bound[0])*.5f);
  float h = MIN(20, (bound[3]-bound[1])*.5f);

  // transform to center of gaussian and round to pixels:
  bound[0] = floorf(i + mu[0]-w);
  bound[1] = floorf(j + mu[1]-h);
  bound[2] = ceilf (i + mu[0]+w);
  bound[3] = ceilf (j + mu[1]+h);

  // normalise
  float weight = 0.0f;
  for(int pj=bound[1];pj<=bound[3];pj++) for(int pi=bound[0];pi<=bound[2];pi++)
  {
    if(pi < 0 || pj < 0 || pi >= wd || pj >= ht) continue;
    float v[2] = {0.5f + pi - i - mu[0], 0.5f + pj - j - mu[1]};
    float vt[2];
    mat2_mulv(B, v, vt);
    float f = det * fakegaussian_pdf(vt[0])*fakegaussian_pdf(vt[1]);
    weight += f;
  }

  weight = 1.0f/weight;
  if(!(weight < FLT_MAX)) return;

  // go through all pixels, multiply by factored covariance matrix, evaluate fake gaussian
  for(int pj=bound[1];pj<=bound[3];pj++) for(int pi=bound[0];pi<=bound[2];pi++)
  {
    if(pi < 0 || pj < 0 || pi >= wd || pj >= ht) continue;
    float v[2] = {0.5f + pi - i - mu[0], 0.5f + pj - j - mu[1]};
    float vt[2];
    mat2_mulv(B, v, vt);
    float f = weight * det * fakegaussian_pdf(vt[0])*fakegaussian_pdf(vt[1]);
    const float col2[3] = {col[0]*f, col[1]*f, col[2]*f};
    filter_box_splat(fb, pi, pj, col2);
  }
}

