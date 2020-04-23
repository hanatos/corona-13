#pragma once

#include "aperture.h"

#include <math.h>
#include <float.h>

// for some jacobian magic, we need to invert a 4x4 matrix:

// fake 5x5 determinant. assumes the wavelength (5th component) stays constant
// and just computes the determinant of the 4x4 subblock.
static inline float lens_det5(const float X[][5])
{
  return
      X[0][3] * X[1][2] * X[2][1] * X[3][0]
    - X[0][2] * X[1][3] * X[2][1] * X[3][0]
    - X[0][3] * X[1][1] * X[2][2] * X[3][0]
    + X[0][1] * X[1][3] * X[2][2] * X[3][0]
    + X[0][2] * X[1][1] * X[2][3] * X[3][0]
    - X[0][1] * X[1][2] * X[2][3] * X[3][0]
    - X[0][3] * X[1][2] * X[2][0] * X[3][1]
    + X[0][2] * X[1][3] * X[2][0] * X[3][1]
    + X[0][3] * X[1][0] * X[2][2] * X[3][1]
    - X[0][0] * X[1][3] * X[2][2] * X[3][1]
    - X[0][2] * X[1][0] * X[2][3] * X[3][1]
    + X[0][0] * X[1][2] * X[2][3] * X[3][1]
    + X[0][3] * X[1][1] * X[2][0] * X[3][2]
    - X[0][1] * X[1][3] * X[2][0] * X[3][2]
    - X[0][3] * X[1][0] * X[2][1] * X[3][2]
    + X[0][0] * X[1][3] * X[2][1] * X[3][2]
    + X[0][1] * X[1][0] * X[2][3] * X[3][2]
    - X[0][0] * X[1][1] * X[2][3] * X[3][2]
    - X[0][2] * X[1][1] * X[2][0] * X[3][3]
    + X[0][1] * X[1][2] * X[2][0] * X[3][3]
    + X[0][2] * X[1][0] * X[2][1] * X[3][3]
    - X[0][0] * X[1][2] * X[2][1] * X[3][3]
    - X[0][1] * X[1][0] * X[2][2] * X[3][3]
    + X[0][0] * X[1][1] * X[2][2] * X[3][3];
}

// fake 5x5 matrix inverse (only inverts the upper left 4x4 sub block).
// last dimension is wavelength which stays unchanged so we leave that part as is.
// this is used by the reverse evaluation newton iteration, to make it match
// the forward evaluation.
static inline void lens_mat5inv(const float X[][5], float R[][5])
{
  const float det = lens_det5(X);
  R[0][0] = ( X[1][2]*X[2][3]*X[3][1] - X[1][3]*X[2][2]*X[3][1] + X[1][3]*X[2][1]*X[3][2] - X[1][1]*X[2][3]*X[3][2] - X[1][2]*X[2][1]*X[3][3] + X[1][1]*X[2][2]*X[3][3] ) / det;
  R[1][0] = ( X[1][3]*X[2][2]*X[3][0] - X[1][2]*X[2][3]*X[3][0] - X[1][3]*X[2][0]*X[3][2] + X[1][0]*X[2][3]*X[3][2] + X[1][2]*X[2][0]*X[3][3] - X[1][0]*X[2][2]*X[3][3] ) / det;
  R[2][0] = ( X[1][1]*X[2][3]*X[3][0] - X[1][3]*X[2][1]*X[3][0] + X[1][3]*X[2][0]*X[3][1] - X[1][0]*X[2][3]*X[3][1] - X[1][1]*X[2][0]*X[3][3] + X[1][0]*X[2][1]*X[3][3] ) / det;
  R[3][0] = ( X[1][2]*X[2][1]*X[3][0] - X[1][1]*X[2][2]*X[3][0] - X[1][2]*X[2][0]*X[3][1] + X[1][0]*X[2][2]*X[3][1] + X[1][1]*X[2][0]*X[3][2] - X[1][0]*X[2][1]*X[3][2] ) / det;

  R[0][1] = ( X[0][3]*X[2][2]*X[3][1] - X[0][2]*X[2][3]*X[3][1] - X[0][3]*X[2][1]*X[3][2] + X[0][1]*X[2][3]*X[3][2] + X[0][2]*X[2][1]*X[3][3] - X[0][1]*X[2][2]*X[3][3] ) / det;
  R[1][1] = ( X[0][2]*X[2][3]*X[3][0] - X[0][3]*X[2][2]*X[3][0] + X[0][3]*X[2][0]*X[3][2] - X[0][0]*X[2][3]*X[3][2] - X[0][2]*X[2][0]*X[3][3] + X[0][0]*X[2][2]*X[3][3] ) / det;
  R[2][1] = ( X[0][3]*X[2][1]*X[3][0] - X[0][1]*X[2][3]*X[3][0] - X[0][3]*X[2][0]*X[3][1] + X[0][0]*X[2][3]*X[3][1] + X[0][1]*X[2][0]*X[3][3] - X[0][0]*X[2][1]*X[3][3] ) / det;
  R[3][1] = ( X[0][1]*X[2][2]*X[3][0] - X[0][2]*X[2][1]*X[3][0] + X[0][2]*X[2][0]*X[3][1] - X[0][0]*X[2][2]*X[3][1] - X[0][1]*X[2][0]*X[3][2] + X[0][0]*X[2][1]*X[3][2] ) / det;

  R[0][2] = ( X[0][2]*X[1][3]*X[3][1] - X[0][3]*X[1][2]*X[3][1] + X[0][3]*X[1][1]*X[3][2] - X[0][1]*X[1][3]*X[3][2] - X[0][2]*X[1][1]*X[3][3] + X[0][1]*X[1][2]*X[3][3] ) / det;
  R[1][2] = ( X[0][3]*X[1][2]*X[3][0] - X[0][2]*X[1][3]*X[3][0] - X[0][3]*X[1][0]*X[3][2] + X[0][0]*X[1][3]*X[3][2] + X[0][2]*X[1][0]*X[3][3] - X[0][0]*X[1][2]*X[3][3] ) / det;
  R[2][2] = ( X[0][1]*X[1][3]*X[3][0] - X[0][3]*X[1][1]*X[3][0] + X[0][3]*X[1][0]*X[3][1] - X[0][0]*X[1][3]*X[3][1] - X[0][1]*X[1][0]*X[3][3] + X[0][0]*X[1][1]*X[3][3] ) / det;
  R[3][2] = ( X[0][2]*X[1][1]*X[3][0] - X[0][1]*X[1][2]*X[3][0] - X[0][2]*X[1][0]*X[3][1] + X[0][0]*X[1][2]*X[3][1] + X[0][1]*X[1][0]*X[3][2] - X[0][0]*X[1][1]*X[3][2] ) / det;

  R[0][3] = ( X[0][3]*X[1][2]*X[2][1] - X[0][2]*X[1][3]*X[2][1] - X[0][3]*X[1][1]*X[2][2] + X[0][1]*X[1][3]*X[2][2] + X[0][2]*X[1][1]*X[2][3] - X[0][1]*X[1][2]*X[2][3] ) / det;
  R[1][3] = ( X[0][2]*X[1][3]*X[2][0] - X[0][3]*X[1][2]*X[2][0] + X[0][3]*X[1][0]*X[2][2] - X[0][0]*X[1][3]*X[2][2] - X[0][2]*X[1][0]*X[2][3] + X[0][0]*X[1][2]*X[2][3] ) / det;
  R[2][3] = ( X[0][3]*X[1][1]*X[2][0] - X[0][1]*X[1][3]*X[2][0] - X[0][3]*X[1][0]*X[2][1] + X[0][0]*X[1][3]*X[2][1] + X[0][1]*X[1][0]*X[2][3] - X[0][0]*X[1][1]*X[2][3] ) / det;
  R[3][3] = ( X[0][1]*X[1][2]*X[2][0] - X[0][2]*X[1][1]*X[2][0] + X[0][2]*X[1][0]*X[2][1] - X[0][0]*X[1][2]*X[2][1] - X[0][1]*X[1][0]*X[2][2] + X[0][0]*X[1][1]*X[2][2] ) / det;
}


// helper function for dumped polynomials to compute integer powers of x:
static inline float lens_ipow(const float x, const int exp)
{
  if(exp == 0) return 1.0f;
  if(exp == 1) return x;
  if(exp == 2) return x*x;
  const float p2 = lens_ipow(x, exp/2);
  if(exp &  1) return x * p2 * p2;
  return p2 * p2;
}

static inline void lens_sphereToCs(const float *inpos, const float *indir, float *outpos, float *outdir, const float sphereCenter, const float sphereRad)
{
  const float normal[3] =
  {
    inpos[0]/sphereRad,
    inpos[1]/sphereRad,
    sqrtf(MAX(0, sphereRad*sphereRad-inpos[0]*inpos[0]-inpos[1]*inpos[1]))/fabsf(sphereRad)
  };
  const float tempDir[3] = {indir[0], indir[1], sqrtf(MAX(0.0, 1.0f-indir[0]*indir[0]-indir[1]*indir[1]))};

  float ex[3] = {normal[2], 0, -normal[0]};
  normalise(ex);
  float ey[3];
  crossproduct(normal, ex, ey);

  outdir[0] = tempDir[0] * ex[0] + tempDir[1] * ey[0] + tempDir[2] * normal[0];
  outdir[1] = tempDir[0] * ex[1] + tempDir[1] * ey[1] + tempDir[2] * normal[1];
  outdir[2] = tempDir[0] * ex[2] + tempDir[1] * ey[2] + tempDir[2] * normal[2];
  outpos[0] = inpos[0];
  outpos[1] = inpos[1];
  outpos[2] = normal[2] * sphereRad + sphereCenter;
}

static inline void lens_csToSphere(const float *inpos, const float *indir, float *outpos, float *outdir, const float sphereCenter, const float sphereRad)
{
  const float normal[3] =
  {
    inpos[0]/sphereRad,
    inpos[1]/sphereRad,
    fabsf((inpos[2]-sphereCenter)/sphereRad)
  };
  float tempDir[3] = {indir[0], indir[1], indir[2]};
  normalise(tempDir);

  float ex[3] = {normal[2], 0, -normal[0]};
  normalise(ex);
  float ey[3];
  crossproduct(normal, ex, ey);
  outdir[0] = dotproduct(tempDir, ex);
  outdir[1] = dotproduct(tempDir, ey);
  outpos[0] = inpos[0];
  outpos[1] = inpos[1];
}

// a few static defines about general lens geometry
#include "init.h"

// evaluates from sensor (in) to outer pupil (out).
// input arrays are 5d [x,y,dx,dy,lambda] where dx and dy are the direction in
// two-plane parametrization (that is the third component of the direction would be 1.0).
// units are millimeters for lengths and micrometers for the wavelength (so visible light is about 0.4--0.7)
// returns the transmittance computed from the polynomial.
static inline float lens_evaluate(const float *in, float *out)
{
  __attribute__ ((unused)) const float x = in[0], y = in[1], dx = in[2], dy = in[3], lambda = in[4];
#include "pt_evaluate.h"
  out[0] = out_x; out[1] = out_y; out[2] = out_dx; out[3] = out_dy;
  return MAX(0.0f, out_transmittance);
}

// evaluates from the sensor (in) to the aperture (out) only
// returns the transmittance.
static inline float lens_evaluate_aperture(const float *in, float *out)
{
  __attribute__ ((unused)) const float x = in[0], y = in[1], dx = in[2], dy = in[3], lambda = in[4];
#include "pt_evaluate_aperture.h"
  out[0] = out_x; out[1] = out_y; out[2] = out_dx; out[3] = out_dy;
  return MAX(0.0f, out_transmittance);
}

// solves for the two directions [dx,dy], keeps the two positions [x,y] and the
// wavelength, such that the path through the lens system will be valid, i.e.
// lens_evaluate_aperture(in, out) will yield the same out given the solved for in.
// in: point on sensor. out: point on aperture.
static inline void lens_pt_sample_aperture(float *in, float *out, float dist)
{
  __attribute__ ((unused)) float out_x = out[0], out_y = out[1], out_dx = out[2], out_dy = out[3], out_transmittance = 1.0f;
  __attribute__ ((unused)) float x = in[0], y = in[1], dx = in[2], dy = in[3], lambda = in[4];
#include "pt_sample_aperture.h"
  // directions may have changed, copy all to be sure.
  out[0] = out_x; out[1] = out_y; out[2] = out_dx; out[3] = out_dy;
  in[0] = x; in[1] = y; in[2] = dx; in[3] = dy;
}

// solves for a sensor position give a scene point and an aperture point
// returns transmittance from sensor to outer pupil
static inline float lens_lt_sample_aperture(
    const float *scene,   // 3d point in scene in camera space
    const float *ap,      // 2d point on aperture (in camera space, z is known)
    float *sensor,        // output point and direction on sensor plane/plane
    float *out,           // output point and direction on outer pupil
    const float lambda)   // wavelength
{
  const float scene_x = scene[0], scene_y = scene[1], scene_z = scene[2];
  const float ap_x = ap[0], ap_y = ap[1];
  float x = 0, y = 0, dx = 0, dy = 0;
#include "lt_sample_aperture.h"
  sensor[0] = x; sensor[1] = y; sensor[2] = dx; sensor[3] = dy; sensor[4] = lambda;
  return MAX(0.0f, out[4]);
}

// jacobian of polynomial mapping sensor to outer pupil. in[]: sensor point/direction/lambda.
static inline void lens_evaluate_jacobian(const float *in, float *J)
{
  __attribute__ ((unused)) const float x = in[0], y = in[1], dx = in[2], dy = in[3], lambda = in[4];
#include "pt_evaluate_jacobian.h"
  J[0]  = dx00; J[1]  = dx01; J[2]  = dx02; J[3]  = dx03; J[4]  = dx04;
  J[5]  = dx10; J[6]  = dx11; J[7]  = dx12; J[8]  = dx13; J[9]  = dx14;
  J[10] = dx20; J[11] = dx21; J[12] = dx22; J[13] = dx23; J[14] = dx24;
  J[15] = dx30; J[16] = dx31; J[17] = dx32; J[18] = dx33; J[19] = dx34;
  J[20] = dx40; J[21] = dx41; J[22] = dx42; J[23] = dx43; J[24] = dx44;
}

static inline float lens_det_sensor_to_outer_pupil(const float *sensor, const float *out, const float focus)
{
  float J[25];
  lens_evaluate_jacobian(sensor, J);
  // only interested in how the direction density at the sensor changes wrt the vertex area on the output
  float T[25] = {
    1., 0., focus, 0., 0.,
    0., 1., 0., focus, 0.,
    0., 0., 1., 0., 0.,
    0., 0., 0., 1., 0.,
    0., 0., 0., 0., 1.};
  float JT[25] = {0.};
  for(int i=2;i<4;i++) // only interested in 2x2 subblock.
    for(int j=0;j<2;j++)
      for(int k=0;k<4;k++)
        JT[i+5*j] += J[k + 5*j] * T[i + 5*k];
  const float det = JT[2] * JT[5+3] - JT[3] * JT[5+2];

  // convert from projected disk to point on hemi-sphere
  const float R = lens_outer_pupil_curvature_radius;
  const float deto = sqrtf(R*R-out[0]*out[0]-out[1]*out[1])/R;
  // there are two spatial components which need conversion to dm:
  const float dm2mm = 100.0f;
  return fabsf(det * deto) / (dm2mm*dm2mm);
}

static inline void lens_evaluate_aperture_jacobian(const float *in, float *J)
{
  __attribute__ ((unused)) const float x = in[0], y = in[1], dx = in[2], dy = in[3], lambda = in[4];
#include "pt_evaluate_aperture_jacobian.h"
  J[0]  = dx00; J[1]  = dx01; J[2]  = dx02; J[3]  = dx03; J[4]  = dx04;
  J[5]  = dx10; J[6]  = dx11; J[7]  = dx12; J[8]  = dx13; J[9]  = dx14;
  J[10] = dx20; J[11] = dx21; J[12] = dx22; J[13] = dx23; J[14] = dx24;
  J[15] = dx30; J[16] = dx31; J[17] = dx32; J[18] = dx33; J[19] = dx34;
  J[20] = dx40; J[21] = dx41; J[22] = dx42; J[23] = dx43; J[24] = dx44;
}

// compute 2x2 determinant of the jacobian mapping position on the aperture to position on the outer pupil.
static inline float lens_det_aperture_to_outer_pupil(const float *sensor, const float *out)
{
  // first get full jacobians from sensor to outer pupil and aperture:
  float Jo[25] = {0.0f}, Ja[25] = {0.0f};
  lens_evaluate_jacobian(sensor, Jo);
  lens_evaluate_aperture_jacobian(sensor, Ja);
  // invert aperture jacobian to map from aperture to sensor:
  float Ja_inv[25];
  lens_mat5inv((const float (*)[5])Ja, (float (*)[5])Ja_inv);
  // multiply to get final mapping:
  float Ja_o[25] = {0.0f};
  for(int i=0;i<2;i++) // only interested in spatial/spatial 2x2 subblock.
    for(int j=0;j<2;j++)
      for(int k=0;k<4;k++)
        Ja_o[i+5*j] += Jo[k + 5*j] * Ja_inv[i + 5*k];
  const float det_ao = Ja_o[0] * Ja_o[5+1] - Ja_o[1] * Ja_o[5+0];
  const float R = lens_outer_pupil_curvature_radius;
  const float deto = R/sqrtf(R*R-out[0]*out[0]-out[1]*out[1]);
  // return determinant:
  return fabsf(deto * det_ao);
}

static inline float lens_det_aperture_to_sensor(const float *sensor, const float focus)
{
  float J[25];
  lens_evaluate_aperture_jacobian(sensor, J);
  // only interested in how the directional density at the sensor changes wrt the vertex area (spatial) at the aperture
  float T[25] = {
    1., 0., focus, 0., 0.,
    0., 1., 0., focus, 0.,
    0., 0., 1., 0., 0.,
    0., 0., 0., 1., 0.,
    0., 0., 0., 0., 1.};
  float JT[25] = {0.};
  for(int i=2;i<4;i++) // only interested in 2x2 subblock.
    for(int j=0;j<2;j++)
      for(int k=0;k<4;k++)
        JT[i+5*j] += J[k + 5*j] * T[i + 5*k];
  const float det = fabsf(JT[2] * JT[5+3] - JT[3] * JT[5+2]);
  // there are two spatial components which need conversion to dm:
  const float dm2mm = 100.0f;
  return dm2mm*dm2mm/det;
}

