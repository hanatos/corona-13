#pragma once
#include <complex.h>
#include "sampler_common.h"
#include "points.h"

#define komplex complex

#ifndef MICRO_MAX_BOUNCES
#define MICRO_MAX_BOUNCES 3
#endif
// #define MICRO_MAX_BOUNCES 1 // only direct regular bsdf
// XXX all but 1 seems utterly broken!
#define MICRO_WALK_MODE 1 // 1: single walk, 2: more grazing walk, 3: bidir + MIS

// this is a parametric header, you can change the behaviour by defining:
// 1) materials (governing fresnel and transmission or extinction)
// MICRO_MATERIAL_DIELECTRIC // dielectric phase function
// MICRO_MATERIAL_CONDUCTOR  // conductors
// MICRO_MATERIAL_DIFFUSE    // diffuse
// 2) slope distributions
// MICRO_SLOPE_GGX
// MICRO_SLOPE_BECKMANN // XXX unimplemented
// 3) height distributions
// MICRO_HEIGHT_UNIFORM
// MICRO_HEIGHT_GAUSSIAN // XXX unimplemented


// return projected roughness for given incident direction and anisotropic roughnesses alpha x,y
// does not depend on sign of wi.
static inline float micro_projected_roughness(const float *wi, const float alpha_x, const float alpha_y)
{
  const float inv_sin_theta2 = 1.0f / (1.0f - wi[2]*wi[2]);
  const float cos_phi2 = wi[0]*wi[0]*inv_sin_theta2;
  const float sin_phi2 = wi[1]*wi[1]*inv_sin_theta2;
  return sqrtf(cos_phi2*alpha_x*alpha_x + sin_phi2*alpha_y*alpha_y);
}

//===================================================================
// microslope distribution:
//===================================================================
#ifdef MICRO_SLOPE_GGX // ggx
static inline float micro_slope_projected_area(const float *wi, const float alpha_x, const float alpha_y)
{
  if(wi[2] >  0.9999f) return 1.0f;
  if(wi[2] < -0.9999f) return 0.0f;

  const float sin2_theta_i = fmaxf(0.0f, 1.0f-wi[2]*wi[2]);
  const float alphai = micro_projected_roughness(wi, alpha_x, alpha_y);
  return 0.5f * (wi[2] + sqrtf(wi[2]*wi[2] + sin2_theta_i*alphai*alphai));
}

static inline float micro_slope_lambda(const float slope, const float roughness)
{
  if(slope >=  1e20f) return 0.0f;
  if(slope <= -1e20f) return 0.0f;
  const float ai = roughness/slope;
  // fprintf(stderr, "our a : %g\n", 1.0/ai);
  // fprintf(stderr, "our lambda: %g\n", 
  //   0.5f*(-1.0f + copysignf(sqrtf(1 + ai*ai), ai)));
  return 0.5f*(-1.0f + copysignf(sqrtf(1 + ai*ai), ai));
}

static inline void micro_slope_sample_11(
    // input
    const float theta_i, // theta incident (normal, wi)
    float U1, float U2, // random numbers
    // output
    float *slope_x,
    float *slope_y)
{
  // special case (normal incidence)
  // XXX FIXME: this threshold is unwisely chosen, it causes inconsistency between sampling and bsdf for normal incidence!
  if(theta_i < 0.0001f)
  {
    const float r = sqrtf(U1/MAX(1e-8f, 1-U1));
    const float phi = 2.0f * M_PI * U2;
    *slope_x = r * cosf(phi);
    *slope_y = r * sinf(phi);
    // fprintf(stderr, "our slope considered degenerate theta %g\n", theta_i);
    return;
  }
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i/cos_theta_i;

#if 0
  // precomputations
  const float a = 1.0f / tan_theta_i;
  const float G1 = 2.0f / (1.0f + sqrtf(1.0f+1.0f/(a*a)));

  // sample slope_x
  const float A = 2.0*U1/G1 - 1.0;
  const float tmp = 1.0 / (A*A-1.0);
  const float B = tan_theta_i;
  const float D = sqrtf(fmaxf(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
  float slope_x_1 = B*tmp - D;
  float slope_x_2 = B*tmp + D;
  // fix cases where B = inf and tmp = 0 and tmp = inf etc.
  if(!(slope_x_1 == slope_x_1)) slope_x_1 = 0.0f;
  if(!(slope_x_2 == slope_x_2)) slope_x_2 = 0.0f;
  assert(slope_x_1 == slope_x_1);
  assert(slope_x_2 == slope_x_2);
  *slope_x = (A < 0 || slope_x_2 > 1.0/tan_theta_i) ? slope_x_1 : slope_x_2;
#endif
	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f*U1/cos_theta_i/c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A*A-1.0f);

	const float D = sqrtf(MAX(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const float slope_x_1 = B*tmp - D;
	const float slope_x_2 = B*tmp + D;
	*slope_x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;
	if(!(*slope_x == *slope_x))
  {
    *slope_x = *slope_y = 0.0f;
    return;
  }
  // fprintf(stderr, "our slope x %g\n", *slope_x);
  // sample slope_y
  float S;
  if(U2 > 0.5)
  {
    S = 1.0;
    U2 = 2.0*(U2-0.5);
  }
  else
  {
    S = -1.0;
    U2 = 2.0*(0.5-U2);
  }
  // improved fit from mitsuba:
  const float z = (U2 * (U2 * (U2 * (-0.365728915865723f) + 0.790235037209296f) -
        0.424965825137544f) + 0.000152998850436920f) /
    (U2 * (U2 * (U2 * (U2 * 0.169507819808272f - 0.397203533833404f) -
                 0.232500544458471f) + 1.0f) - 0.539825872510702f);
  // fprintf(stderr, "our z %g\n", z);
  *slope_y = S * z * sqrtf(1.0+*slope_x * *slope_x);
  assert(*slope_y == *slope_y);
}

static inline float micro_slope_P22(const float slope_x, const float slope_y, const float alpha_x, const float alpha_y)
{
  const float tmp = 1.0f + slope_x*slope_x/(alpha_x*alpha_x) + slope_y*slope_y/(alpha_y*alpha_y);
  return 1.0f / (M_PI * alpha_x * alpha_y) / (tmp * tmp);
}
#endif
#ifdef MICRO_SLOPE_BECKMANN // beckmann
#endif


//===================================================================
// microsurface height distribution:
//===================================================================
#if 1 // uniform
// uniform height distribution:
static inline float micro_height_C1(const float h)
{
  return fminf(1.0, fmaxf(0.0, 0.5f*(h+1.0f)));
}
static inline float micro_height_inv_C1(const float u)
{
  return fmaxf(-1.0f, fminf(1.0, 2.0f*u-1.0f));
}
#endif


//===================================================================
// microfacet slope distribution evaluation and sampling
//===================================================================

// evaluate slope distribution D at half vector h in tangent space
static inline float micro_slope_D(const float *h, const float alpha_x, const float alpha_y)
{
  if(h[2] <= 0.0) return 0.0f;
  const float slope_x = -h[0]/h[2];
  const float slope_y = -h[1]/h[2];
  return micro_slope_P22(slope_x, slope_y, alpha_x, alpha_y) / (h[2]*h[2]*h[2]*h[2]);
}

// evaluate slope distribution of visible normals for tangent space half vector h
// and incident direction wi pointing away from the surface point
static inline float micro_slope_D_wi(const float *wi, const float *h, const float alpha_x, const float alpha_y)
{
  if(h[2] <= 0.0) return 0.0f;
  const float projected_area = micro_slope_projected_area(wi, alpha_x, alpha_y);
  if(projected_area == 0) return 0;
  const float c = 1.0f / projected_area;
  return c * fmaxf(0.0f, dotproduct(wi, h)) * micro_slope_D(h, alpha_x, alpha_y);
}

static inline void micro_slope_sample_D(const float *wi, const float r1, const float r2, const float alpha_x, const float alpha_y, float *h)
{
  // fprintf(stderr, "our sample D input %g %g %g \n", wi[0], wi[1], wi[2]);
  // 1. stretch omega_i
  float wi_[3] = { alpha_x * wi[0], alpha_y* wi[1], wi[2] };
  normalise(wi_);
  // get polar coordinates of omega_i_
  float theta = 0.0f;
  float sin_phi = 0.0f, cos_phi = 1.0f;
  if (wi_[2] < 0.99999)
  {
    const float len = sqrtf(wi_[0]*wi_[0] + wi_[1]*wi_[1]);
    // fprintf(stderr, "our theta_i %g\n", acosf(wi_[2]));
    theta = acosf(wi_[2]);//len/wi_[2];
    sin_phi = wi_[1]/len;
    cos_phi = wi_[0]/len;
  }
  // 2. sample P22_{wi}(x_slope, y_slope, 1, 1)
  float slope_x, slope_y;
  micro_slope_sample_11(theta, r1, r2, &slope_x, &slope_y);
  // 3. rotate
#if 1
  float tmp = cos_phi*slope_x - sin_phi*slope_y;
  slope_y = sin_phi*slope_x + cos_phi*slope_y;
  slope_x = tmp;
#else
  const float phi = atan2f(wi_[1], wi_[0]);
  float tmp = cosf(phi)*slope_x - sinf(phi)*slope_y;
  slope_y = sinf(phi)*slope_x + cosf(phi)*slope_y;
  slope_x = tmp;
#endif
  // 4. unstretch
  slope_x *= alpha_x;
  slope_y *= alpha_y;
  // fprintf(stderr, "our unstretched slope %g %g\n", slope_x, slope_y);
  // 5. compute normal
  float inv_h = 1.0f/sqrtf(slope_x*slope_x + slope_y*slope_y + 1.0f);
  h[0] = -slope_x*inv_h;
  h[1] = -slope_y*inv_h;
  h[2] = inv_h;
  if(!(inv_h > 0.0))
  {
    h[0] = h[2] = 0.0f;
    h[1] = 1.0f;
  }
}

// fwd declare correlated heights G2 and G1 averaged over height
static inline float micro_G2(const float *wi, const float *wo, const float alpha_x, const float alpha_y);
static inline float micro_G1_avg(const float slope, const float roughness);


//===================================================================
// surface material:
//===================================================================
#ifdef MICRO_MATERIAL_CONDUCTOR // conductors:
static inline float micro_surface_fresnel(const float n1, const komplex float n2, const float cosr)
{
  const komplex float cost = csqrtf(1.0f - (n1/n2)*(n1/n2) * (1.0f - cosr*cosr));
  // fresnel for unpolarized light:
  const float Rs = cabsf((n1*cosr - n2*cost)/(n1*cosr + n2*cost));
  const float Rp = cabsf((n1*cost - n2*cosr)/(n1*cost + n2*cosr));
  return fminf(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}
// returns new outside flag
static inline float micro_sample_phase_function(
    float *dir,          // ray direction (pointing towards intersection), will be overwritten
    const float r0,      // 3 random numbers: direction and layer selection
    const float r1,
    const float r2,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float n1,            // ior of the R material
    komplex float n2,    // ior of the T material
    int *inside)         // flag whether inside the material
{
  float h[3], wi[3];
  // wi for sampling needs to point away from surface!
  for(int k=0;k<3;k++) wi[k] = -dir[k];
  micro_slope_sample_D(wi, r1, r2, alpha_x, alpha_y, h);
  const float dot = dotproduct(dir, h);
  for(int k=0;k<3;k++) dir[k] = dir[k] - 2.0f*h[k] * dot;
  return micro_surface_fresnel(n1, n2, fabsf(dot));
}

static inline float micro_eval_phase_function(
    const float *wi,     // incoming direction, pointing towards surface
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    const float n1,
    const komplex float n2,
    const int inside)    // inside means below the tangent-space surface at n=(0,0,1)
{
  float h[3];
  for(int k=0;k<3;k++) h[k] = - wi[k] + wo[k];
  normalise(h);
  const float cosr = fabsf(dotproduct(wi, h));
  const float F = micro_surface_fresnel(n1, n2, cosr);

  const float nwi[3] = {-wi[0], -wi[1], -wi[2]};
  const float D_wi = micro_slope_D_wi(nwi, h, alpha_x, alpha_y);
  return F * D_wi / (4.0f * cosr);
}

// same as eval phase func, but without fresnel for conductors.
static inline float micro_pdf_phase_function(
    const float *wi,     // ray direction (pointing towards intersection)
    const float *wo,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float n1,
    komplex float n2,    // ior of T material
    int inside)          // flag whether inside the material
{
  float h[3];
  for(int k=0;k<3;k++) h[k] = - wi[k] + wo[k];
  normalise(h);
  const float cosr = fabsf(dotproduct(wi, h));
  const float nwi[3] = {-wi[0], -wi[1], -wi[2]};
  const float D_wi = micro_slope_D_wi(nwi, h, alpha_x, alpha_y);
  return D_wi / (4.0f * cosr);
}

static inline float micro_eval(
    const float *wi,     // incoming direction, pointing towards surface, wi[2] < 0
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    float n1,
    komplex float n2)    // ior of T material
{
  float h[3];
  for(int k=0;k<3;k++) h[k] = - wi[k] + wo[k];
  normalise(h);
  const float cosr = fabsf(dotproduct(wi, h));
  const float G2 = micro_G2(wi, wo, alpha_x, alpha_y);
  const float F = micro_surface_fresnel(n1, n2, cosr);
  const float D = micro_slope_D(h, alpha_x, alpha_y);
  return F * D * G2 / (4.0f * cosr);
}
#endif
#ifdef MICRO_MATERIAL_DIELECTRIC // dielectric phase function:
static inline float micro_surface_cost(const float n1, const float n2, const float cosr)
{
  const float cost2 = 1.0f - (n1/n2)*(n1/n2) * (1.0f - cosr*cosr);
  return (cost2 <= 0.0f) ? 0.0f : sqrtf(cost2);
}
static inline float micro_surface_fresnel(const float n1, const float n2, const float cosr, const float cost)
{
  if(cost <= 0.0f) return 1.0f; // total inner reflection
  // fresnel for unpolarized light:
  const float Rs = (n1*cosr - n2*cost)/(n1*cosr + n2*cost);
  const float Rp = (n1*cost - n2*cosr)/(n1*cost + n2*cosr);
  return fminf(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}
// returns new outside flag
static inline float micro_sample_phase_function(
    float *dir,          // ray direction (pointing towards intersection), will be overwritten
    const float r0,      // 3 random numbers: direction and layer selection
    const float r1,
    const float r2,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float n1,
    komplex float n2,    // ior of T material
    int *inside)         // flag whether inside the material
{
  // fprintf(stderr, "ours got n1 n2 %g %g\n", n1, n2);
  const float eta = *inside ? n1/crealf(n2) : crealf(n2)/n1;
  // wi for sampling needs to point away from surface!
  float h[3], wi[3] = {-dir[0], -dir[1], -dir[2]};
  // if(*inside)
  if(dir[2] > 0)
  {wi[0] = -wi[0], wi[1] = -wi[1]; wi[2] = -wi[2];}
  // if(*inside) wi[2] = -wi[2];
  micro_slope_sample_D(wi, r1, r2, alpha_x, alpha_y, h);
  if(dir[2] > 0)
  // if(*inside)
  {h[0] = -h[0]; h[1] = -h[1]; h[2] = -h[2];};
  // if(*inside) h[2] = -h[2];
  const float cosr = -dotproduct(h, dir);
  // fprintf(stderr, "our outside %d\n", !*inside);
  // fprintf(stderr, "our cosr %g\n", cosr);
  // fprintf(stderr, "our h and eta %g %g %g  %g\n", h[0], h[1], h[2], eta);
  // if(cosr < 1e-9f) return 0.0f;
  const float cost = micro_surface_cost(1.0, eta, cosr);
  const float F = micro_surface_fresnel(1.0, eta, cosr, cost);
  // fprintf(stderr, "our F %g cost %g\n", F, cost);
  if(r0 < F)
  { // reflect
    // fprintf(stderr, "ours reflects\n");
    for(int k=0;k<3;k++) dir[k] = dir[k] + 2.0f*h[k] * cosr;
    // fprintf(stderr, "ours reflects %g %g %g\n", dir[0], dir[1], dir[2]);
  }
  else
  { // transmit
    // fprintf(stderr, "ours transmits\n");
    // const float f = eta*cosr - cost;
    // for(int k=0;k<3;k++) dir[k] = dir[k]*eta + f * h[k];
    // for(int k=0;k<3;k++) dir[k] = h[k] * (cosr/eta - cost) - wi[k]/eta;
    float f = cosr/eta - cost;
    for(int k=0;k<3;k++) dir[k] = h[k] * f + dir[k]/eta;

    normalise(dir);
    *inside = !*inside;
  }
  // dielectric does not dissipate energy:
  return 1.0f;
}

static inline float micro_eval_phase_function(
    const float *wi,     // incoming direction, pointing towards surface
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    float n1,
    komplex float n2,    // ior of T material
    const int inside)    // inside means below the tangent-space surface at n=(0,0,1)
{
  int reflect = (wo[2] > 0 && !inside) || (wo[2] < 0 && inside);
  // int reflect = (wo[2] > 0) == !inside;
  float eta = inside ? crealf(n2)/n1 : n1/crealf(n2);
  float h[3];
  if(reflect)
  {
    for(int k=0;k<3;k++) h[k] = - wi[k] + wo[k];
  }
  else
  {
    for(int k=0;k<3;k++) h[k] = eta * wi[k] - wo[k];
    // h will point into the thinner medium, but we'd like it to point into the macroscopic R hemisphere
    if(eta > 1.0) for(int k=0;k<3;k++) h[k] = -h[k]; 
  }

  normalise(h);
  if(!(h[0] == h[0])) return 0.0f;
  const float cosr = -dotproduct(wi, h);
  if(cosr < 1e-4f) return 0.0f;
  const float cost = micro_surface_cost(eta, 1.0, cosr);
  if(cost == 0.0 && !reflect) return 0.0f;
  const float F = micro_surface_fresnel(eta, 1.0, cosr, cost);
  const float nwi[3] = {-wi[0], -wi[1], -wi[2]};
  const float D_wi = micro_slope_D_wi(nwi, h, alpha_x, alpha_y);
  if(reflect)
  {
    // if(wo[2] < 0.0) return 0.0f;
    assert(0.25f * D_wi / cosr * F >= 0.0);
    return 0.25f * D_wi / cosr * F;
  }

  const float costt = dotproduct(wo, h);
  if(costt >= 0) return 0.0f;
  // if(wo[2] > 0.0) return 0.0f;

  // jacobian half vector to outgoing solid angle
  const float den = eta * cosr - cost;
  const float dh_dw = cost/(den*den);
  assert((1.0f-F) * D_wi * dh_dw >= 0.0);
  return (1.0f-F) * D_wi * dh_dw;
}

static inline float micro_pdf_phase_function(
    const float *wi,     // ray direction (pointing towards intersection)
    const float *wo,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float n1,
    komplex float n2,    // ior of T material
    int inside)          // flag whether inside the material
{
  return micro_eval_phase_function(wi, wo, alpha_x, alpha_y, n1, n2, inside);
}

// returns bsdf of first bounce * cos(wo)
static inline float micro_eval(
    const float *wi,     // incoming direction, pointing towards surface, wi[2] < 0
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    float n1,
    komplex float n2)    // ior of T material
{
  int reflect = wo[2] > 0;
  float eta = n1/crealf(n2);
  float h[3];
  if(reflect)
    for(int k=0;k<3;k++) h[k] = - wi[k] + wo[k];
  else
  {
    for(int k=0;k<3;k++) h[k] = eta * wi[k] - wo[k];
    // h will point into the thinner medium, but we'd like it to point into the R hemisphere with wi
    if(eta > 1.0) for(int k=0;k<3;k++) h[k] = -h[k]; 
  }

  normalise(h);
  const float cosr = -dotproduct(wi, h);
  if(cosr < 1e-4f) return 0.0f;
  const float cost = micro_surface_cost(eta, 1.0, cosr);
  if(cost == 0.0f && !reflect) return 0.0f;
  const float G2 = micro_G2(wi, wo, alpha_x, alpha_y);
  const float F = micro_surface_fresnel(eta, 1.0, cosr, cost);
  const float D = micro_slope_D(h, alpha_x, alpha_y);
  if(reflect)
    return F * D * G2 / (4.0f * cosr);

  const float ct = dotproduct(h, wo);
  if(ct >= 0) return 0.0f;
  // jacobian half vector to outgoing solid angle
  const float den = eta * cosr - cost;
  const float dh_dw = cost/(den*den);
  assert( (1.0f-F) * D * G2 * dh_dw >= 0.0);
  return (1.0f-F) * D * G2 * dh_dw * fabsf(cosr/wi[2]);
}
#endif
#ifdef MICRO_MATERIAL_DIFFUSE // diffuse micro facets
static inline float micro_sample_phase_function(
    float *dir,          // ray direction (pointing towards intersection), will be overwritten
    const float r0,      // 3 random numbers: direction and layer selection
    const float r1,
    const float r2,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float albedo,        // space for n1 abused for the diffuse albedo
    komplex float n2,    // ior of T material
    int *inside)         // flag whether inside the material
{
  float wom[3], a[3], b[3], h[3], wi[3] = {-dir[0], -dir[1], -dir[2]};
  // wi for sampling needs to point away from surface!
  micro_slope_sample_D(wi, r1, r2, alpha_x, alpha_y, h);
  const float cosr = -dotproduct(h, dir);
  if(cosr < 1e-9f) return 0.0f;
  // diffuse reflection:
  get_onb(h, a, b);
  sample_cos(wom, wom+1, wom+2, r1, r2);
  for(int k=0;k<3;k++) dir[k] = wom[0]*a[k] + wom[1]*b[k] + wom[2]*h[k];
  return albedo;
}

static inline float micro_eval_phase_function(
    const float *wi,     // incoming direction, pointing towards surface
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    float albedo,
    komplex float n2,    // ior of T material
    const int inside)    // inside means below the tangent-space surface at n=(0,0,1)
{
  const float r1 = points_rand(rt.points, common_get_threadid());
  const float r2 = points_rand(rt.points, common_get_threadid());
  float h[3], nwi[3] = {-wi[0], -wi[1], -wi[2]};
  micro_slope_sample_D(nwi, r1, r2, alpha_x, alpha_y, h);
  return 1.0/M_PI * albedo * fmaxf(0.0, -dotproduct(wi, h));
}

static inline float micro_pdf_phase_function(
    const float *wi,     // ray direction (pointing towards intersection)
    const float *wo,
    const float alpha_x, // roughnesses
    const float alpha_y,
    float albedo,
    komplex float n2,    // ior of T material
    int inside)          // flag whether inside the material
{
  return micro_eval_phase_function(wi, wo, alpha_x, alpha_y, albedo, n2, inside);
}

static inline float micro_eval(
    const float *wi,     // incoming direction, pointing towards surface, wi[2] < 0
    const float *wo,     // outgoing direction, pointing away from surface
    const float alpha_x, // anisotropic roughnesses
    const float alpha_y,
    float albedo,
    komplex float n2)    // ior of T material
{
  const float r1 = points_rand(rt.points, common_get_threadid());
  const float r2 = points_rand(rt.points, common_get_threadid());
  float h[3], nwi[3] = {-wi[0], -wi[1], -wi[2]};
  micro_slope_sample_D(nwi, r1, r2, alpha_x, alpha_y, h);
  // shadowing given masking:
  const float slope_i = fabsf(wi[2]) / sqrtf(fmaxf(0.0, 1.0 - wi[2]*wi[2]));
  const float roughness_i = micro_projected_roughness(wi, alpha_x, alpha_y);
  const float G2_given_G1 = micro_G2(wi, wo, alpha_x, alpha_y)/micro_G1_avg(slope_i, roughness_i);
  return 1.0/M_PI * albedo * fmaxf(0.0, -dotproduct(wi, h)) * G2_given_G1;
}
#endif


//===================================================================
// height and distance sampling
//===================================================================

static inline float micro_G1(const float slope, const float roughness, const float h0)
{
  if(slope > 1e20f)  return 1.0f;
  if(!(slope > 0.0)) return 0.0f;
  return powf(micro_height_C1(h0), micro_slope_lambda(slope, roughness));
}
static inline float micro_G1_avg(const float slope, const float roughness)
{
  if(slope > 1e20f) return 1.0f;
  if(slope <= 0.0) return 0.0f;
  return 1.0f/(1.0f + micro_slope_lambda(slope, roughness));
}
static inline float micro_mu_t(const float slope, const float roughness, const float h0, const float tau)
{
  return micro_slope_lambda(slope, roughness) * logf(micro_height_C1(h0 + slope*tau));
}
static inline float micro_inv_mu_t(const float slope, const float roughness, const float h0, const float V)
{
  return 1.0/slope * (micro_height_inv_C1(expf(V/micro_slope_lambda(slope, roughness))) - h0);
}

static inline double micro_abgam(double x)
{
  const double gam0 = 1./ 12.;
  const double gam1 = 1./ 30.;
  const double gam2 = 53./ 210.;
  const double gam3 = 195./ 371.;
  const double gam4 = 22999./ 22737.;
  const double gam5 = 29944523./ 19733142.;
  const double gam6 = 109535241009./ 48264275462.;
  return 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
    + gam0/(x + gam1/(x + gam2/(x + gam3/(x + gam4 /
              (x + gam5/(x + gam6/x))))));
}

static inline double micro_gamma(double x)
{
  return exp (micro_abgam (x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
}

static inline double micro_beta(double m, double n)
{
  return (micro_gamma (m)*micro_gamma (n)/micro_gamma (m + n));
}

// G term with height correlation. costly for T. might want to go with
// teleportation and decorrelated height for G1_i and G1_o instead.
static inline float micro_G2(
    const float *wi,           // incident direction in tangent space, pointing towards surface   (wi[2] < 0)
    const float *wo,           // outgoing direction in tangent space, pointing away from surface (wo[2] > 0 => R wo[2] < 0 => T)
    const float alpha_x,       // roughnesses
    const float alpha_y)
{
  if(wi[2] >= 0) return 0.0f; // should not be used this way
  if(fabsf(wo[2]) < 1e-3f) return 0.0f; // all masked
  const float slope_i = fabsf(wi[2]) / sqrtf(fmaxf(0.0, 1.0 - wi[2]*wi[2]));
  const float slope_o = fabsf(wo[2]) / sqrtf(fmaxf(0.0, 1.0 - wo[2]*wo[2]));
  const float roughness_i = micro_projected_roughness(wi, alpha_x, alpha_y);
  const float roughness_o = micro_projected_roughness(wo, alpha_x, alpha_y);
  if(fabsf(fabsf(wi[2])-1.0f) < 1e-4f)
  {
    if(fabsf(fabsf(wo[2])-1.0f) < 1e-4f) return 1.0f;
    return micro_G1_avg(slope_o, roughness_o);
  }
  else if(fabsf(fabsf(wo[2])-1.0f) < 1e-4f)
    return micro_G1_avg(slope_i, roughness_i);

  const float l_i = micro_slope_lambda(slope_i, roughness_i), l_o = micro_slope_lambda(slope_o, roughness_o);
  if(wo[2] > 0) // reflection:
    return 1.0/(1.0f + l_i + l_o);
  else // transmission
  {
    const float G2 = micro_beta(l_i+1, l_o+1);
    if(!(G2 >= 0)) return 0.0f; // numerical satan got us.
    assert(G2 >= 0);
    assert(G2 <= 1);
    return G2;
  }
    // or with libm: return tgammaf(l_i+1) * tgammaf(l_o+1) / tgammaf(l_i+l_o+2);
}

// sample new height for given direction wo (tangent space) and random number
static inline float micro_sample_height(const float *wo_, float h0, float alpha_x, float alpha_y, float rand, int inside)
{
  float wo[3] = {
    inside ? -wo_[0] : wo_[0], 
    inside ? -wo_[1] : wo_[1], 
    inside ? -wo_[2] : wo_[2]};
  if(inside) h0 = -h0;
  if(wo[2] >  0.9999f) return FLT_MAX;
  if(wo[2] < -0.9999f) return (inside ? -1 : 1) * micro_height_inv_C1(rand * micro_height_C1(h0));
  if(fabsf(wo[2]) < 1e-4f) return (inside ? -1 : 1) * h0;

  const float slope = wo[2] / sqrtf(fmaxf(0.0, 1.0 - wo[2]*wo[2]));

  const float roughness = micro_projected_roughness(wo, alpha_x, alpha_y);
  if(rand > 1.0f - micro_G1(slope, roughness, h0))
    return FLT_MAX;
  return (inside ? -1 : 1) *
    micro_height_inv_C1(
        micro_height_C1(h0) / powf(1.0f-rand, 1.0f/
        micro_slope_lambda(slope, roughness)));
}

static inline float micro_sample_height2(const float *wo, float *h1, float alpha_x, float alpha_y, float rand, int inside)
{
  float h0 = *h1;
  if(inside) h0 = -h0;
  if(wo[2] >  0.9999f) { *h1 = FLT_MAX; return 1.0; }
  if(wo[2] < -0.9999f) { *h1 = (inside ? -1 : 1) * micro_height_inv_C1(rand * micro_height_C1(h0)); return 1.0;}
  if(fabsf(wo[2]) < 1e-4f) { *h1 = (inside ? -1 : 1) * h0; return 1.0; }

  // if inside the material, flip signs of slope and h0
  const float slope = (inside ? -1.0 : 1.0) * wo[2] / sqrtf(fmaxf(0.0, 1.0 - wo[2]*wo[2]));

  const float roughness = micro_projected_roughness(wo, alpha_x, alpha_y);
  const float G1 = micro_G1(slope, roughness, h0);
  const float mu0 = micro_mu_t(slope, roughness, h0, 0);
  const float tau = micro_inv_mu_t(slope, roughness, h0, mu0 - logf(1.0f-rand * (1.0-G1)));
  *h1 = inside ? -h0 - slope * tau : h0 + slope * tau;
  return 1.0-G1;
}



//===================================================================
// combined sampling and evaluation
//===================================================================

// tiny encryption algorithm:
static inline void
micro_sample_two_float(unsigned int v0, unsigned int v1, float *out)
{
  const unsigned int key0 = 0xa341316c;
  const unsigned int key1 = 0xc8013ea4;
  const unsigned int key2 = 0xad90777d;
  const unsigned int key3 = 0x7e95761e;
  unsigned int sum = 0;
  unsigned int delta = 0x9e3779b9;

  // this can be lowered for performance,
  // 16 results in reasonably good points.
  const unsigned int rounds = 16;
  for(unsigned int i = 0; i < rounds; i++)
  {
    sum += delta;
    v0 += ((v1 << 4) + key0) ^ (v1 + sum) ^ ((v1 >> 5) + key1);
    v1 += ((v0 << 4) + key2) ^ (v0 + sum) ^ ((v0 >> 5) + key3);
  }
  // convert to float in [0, 1)
  uint32_t v00 = 0x3f800000 | (v0>>9);
  uint32_t v11 = 0x3f800000 | (v1>>9);
  out[0] = (*(float*)&v00) - 1.0f;
  out[1] = (*(float*)&v11) - 1.0f;
}


static inline float micro_multiple_eval_walk(
    const float *wi,                 // incident direction pointing towards surface, in tangent space
    const float *wo,                 // outgoing direction pointing away from surface, in tangent space
    const float alpha_x,             // anisotropic roughness
    const float alpha_y,
    const float n1,
    const komplex float n2,          // ior with extinction
    int mis,                         // whether to use an mis weight (needs two walks)
    int hash)                        // random seed
{
  float height = 1.0f + micro_height_inv_C1(0.999f);

  float sum = 0.0f;

  // int hash = (int)(points_rand(rt.points, common_get_threadid())*4000000);
  float r0[2], r1[2];
  int seed = 0;
  int inside = 0;
  float throughput = 1.0f;
  float dir[3];
  for(int k=0;k<3;k++) dir[k] = wi[k];

  // slope for next event estimation:
  const float slope_o = wo[2] / sqrtf(fmaxf(0.0, 1.0 - wo[2]*wo[2]));
  const float roughness_o = micro_projected_roughness(wo, alpha_x, alpha_y);

  float pdf_fwd = 1.0f;

  // loop over multiple bounces:
  for(int i=0;i<MICRO_MAX_BOUNCES;i++) // last chance to escape: i = MICRO_MAX_BOUNCES-1
  {
    // draw random numbers for height and R/T
    micro_sample_two_float(hash, seed++, r0);
    // draw random numbers for direction sampling
    micro_sample_two_float(hash, seed++, r1);

    height = micro_sample_height(dir, height, alpha_x, alpha_y, r0[0], inside);
    if(height == FLT_MAX) break;

    if(i)
    { // kill first bounce treated above
      const float G1 = (wo[2] > 0) ? micro_G1(slope_o, roughness_o, height) : micro_G1(-slope_o, roughness_o, -height);
      const float phase_func = micro_eval_phase_function(dir, wo, alpha_x, alpha_y, n1, n2, inside);
      if(mis) // MIS
        sum += throughput * (i ? (pdf_fwd > 0.0f ? pdf_fwd/(pdf_fwd + phase_func) : 0.0) : 0.5) * phase_func * G1;
      else // plain
        sum += throughput * phase_func * G1;
    }

    // overwrite dir to hold new direction, and possibly flip inside bit
    throughput *= micro_sample_phase_function(dir, r0[1], r1[0], r1[1], alpha_x, alpha_y, n1, n2, &inside);
    if(throughput <= 0.0f) break;
    if(mis && i==0) // only compute that first time, so input dir==wi and inside==0
      pdf_fwd = micro_eval_phase_function(wi, dir, alpha_x, alpha_y, n1, n2, 0);
  }
  return sum;
}

static inline float micro_multiple_eval(
    const float *wi,                 // incident direction pointing towards surface, in tangent space
    const float *wo,                 // outgoing direction pointing away from surface, in tangent space
    const float alpha_x,             // anisotropic roughness
    const float alpha_y,
    const float n1,
    const komplex float n2,          // ior with extinction
    const int hash)                  // random seed
{
  if(fabsf(wo[2]) < 1e-9f) return 0.0f;
  // deterministically eval first bounce
  // float first = 0;// XXX indirect only micro_eval(wi, wo, alpha_x, alpha_y, n1, n2)/fabsf(wo[2]);
  float first = micro_eval(wi, wo, alpha_x, alpha_y, n1, n2)/fabsf(wo[2]);
#if MICRO_WALK_MODE == 1
  // variant 1: only one walk
  return first + micro_multiple_eval_walk(wi, wo, alpha_x, alpha_y, n1, n2, 0, hash)/fabsf(wo[2]);
#elif MICRO_WALK_MODE == 2
  // variant 1: choose more grazing walk
  // choose random walk with lower variance:
  // if(wi[2] < 0.0f && <== this should always be the case
  assert(wi[2] < 0.0);
  if(wo[2] > 0.0f)
  { // reflect event
    if(-wi[2] > wo[2])
    { // outgoing direction is more grazing than incoming
      const float wi2[3] = {-wo[0], -wo[1], -wo[2]};
      const float wo2[3] = {-wi[0], -wi[1], -wi[2]};
      return first +
        micro_multiple_eval_walk(wi2, wo2, alpha_x, alpha_y, n1, n2, 0, hash)/fabsf(wo2[2]);
    }
  }
  else if(-wi[2] > -wo[2])
  { // transmit event, outgoing dir is more grazing
    const float wi2[3] = {wo[0], wo[1], wo[2]};
    const float wo2[3] = {wi[0], wi[1], wi[2]};
    return first +
      micro_multiple_eval_walk(wi2, wo2, alpha_x, alpha_y, crealf(n2), n1, 0, hash)/fabsf(wo2[2]);
  }
  return first + micro_multiple_eval_walk(wi, wo, alpha_x, alpha_y, n1, n2, 0, hash)/fabsf(wo[2]);
#elif MICRO_WALK_MODE == 3
  // variant 3: do both walks and weight with MIS
  float sum = micro_multiple_eval_walk(wi, wo, alpha_x, alpha_y, n1, n2, 1, hash)/fabsf(wo[2]);
  const float wi2[3] = {-wo[0], -wo[1], -wo[2]};
  const float wo2[3] = {-wi[0], -wi[1], -wi[2]};
  if(wo[2] > 0.0f) // reflect
    sum += micro_multiple_eval_walk(wi2, wo2, alpha_x, alpha_y, n1, n2, 1, hash)/fabsf(wo2[2]);
  else // transmit
    sum += micro_multiple_eval_walk(wi2, wo2, alpha_x, alpha_y, crealf(n2), n1, 1, hash)/fabsf(wo2[2]);
  return first + sum;
#else
#error "MICRO_WALK_MODE has to be 1, 2, or 3!"
#endif
}

// pdf.
static inline float micro_multiple_pdf(
    const float *wi,                 // incident direction pointing towards surface, in tangent space
    const float *wo,                 // outgoing direction pointing away from surface, in tangent space
    const float alpha_x,             // anisotropic roughness
    const float alpha_y,
    const float n1,
    const komplex float n2,          // ior with extinction
    int hash)                        // random seed
{
  assert(wi[2] <= 0);
  if(fabsf(wo[2]) < 1e-9f) return 0.0f;
  return micro_multiple_eval(wi, wo, alpha_x, alpha_y, n1, n2, hash);
  // deterministic something..:
  // return .5/M_PI + .5f * micro_pdf_phase_function(wi, wo, alpha_x, alpha_y, n1, n2, 0)/fabsf(wo[2]);
}

static inline float micro_multiple_sample(
    const float *wi,         // incident direction pointing towards surface, in tangent space
    float *wo,
    const float alpha_x,     // anisotropic roughness
    const float alpha_y,
    const float n1,
    const komplex float n2,  // ior with extinction
    int hash,                // random seed
    const float rand0,
    const float rand1,
    const float rand2)
{
  float height = 1.0f + micro_height_inv_C1(0.999f);

  float r0[2], r1[2];
  // XXX debug:
  // hash = (int)(points_rand(rt.points, common_get_threadid())*4000000);
  int seed = 0;
  int inside = 0;
  float throughput = 1.0f;
  for(int k=0;k<3;k++) wo[k] = wi[k];

  // loop over multiple bounces:
  for(int i=0;i<=MICRO_MAX_BOUNCES;i++) // last chance to escape: i = MICRO_MAX_BOUNCES-1
  {
    // draw random numbers for height and R/T
    micro_sample_two_float(hash, seed++, r0);
    // draw random numbers for direction sampling
    micro_sample_two_float(hash, seed++, r1);

    height = micro_sample_height(wo, height, alpha_x, alpha_y, r0[0], inside);
    if(height == FLT_MAX)
    {
      assert(i); // should never happen when pointing down
      // check correct sidedness:
      if(wo[2] <= 0 && !inside) return 0.0f;
      if(wo[2] >= 0 &&  inside) return 0.0f;
      // if(i==1) return 0.0; // XXX indirect only
      return throughput;
    }

    // overwrite wo to hold new direction, and possibly flip inside bit
    if(i) throughput *= micro_sample_phase_function(wo, r0[1], r1[0], r1[1], alpha_x, alpha_y, n1, n2, &inside);
    else  throughput *= micro_sample_phase_function(wo, rand0, rand1, rand2, alpha_x, alpha_y, n1, n2, &inside);
    if(throughput <= 0.0f) return 0.0f;
  }
  return 0.0f; // cut off higher orders
  // TODO: i suppose we could replace the last bounce by evaluation of G1
}

