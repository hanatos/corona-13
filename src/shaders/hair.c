/*
    This file is part of corona-13.

    copyright (c) 2015 johannes hanika

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <assert.h>
#include "corona_common.h"
#include "spectrum.h"
#include "shader.h"
#include "prims.h"
#include "pnoise.h"
#include "sampler_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct hair_t
{
  float eumelanin;
  float pheomelanin;
}
hair_t;

// TODO: special case: lambertian fiber (simple case from the paper):
// rd * pdf
// pdf = 1/(4*pi) * ((-phi + pi) * cosf(phi) + sinf(phi))
// and sample diffuse around the cylinder


// input parameters:
// ior,                       n:          1.55
// scales tilt,               alpha_R:    -2..-10 degrees
// internal absorption,       mu_a:       0.2..inf (use {eu,pheo}melanin concentrations)
// variance R,                v_R:        5-10 deg


// n: ior ratio = n_inside / n_outside = n_inside for vacuum,
// cosr: cosine of incoming direction (reflection cosine)
static inline float fresnel(const float n, const float cosr)
{
  const float cost2 = 1.0f - 1.0f/(n*n) * (1.0f - cosr*cosr);
  if(cost2 <= 0.0f) return 1.0f; // total inner reflection
  const float cost = sqrtf(cost2);
  // fresnel for unpolarized light:
  const float Rs = (cosr - n*cost)/(cosr + n*cost);
  const float Rp = (cost - n*cosr)/(cost + n*cosr);
  return fminf(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}

// attenuation due to absorption and fresnel
static inline float A(
    const int p,
    const float h,
    const float n,
    const float np1,
    const float mu_a,
    const float cos_theta_d,
    const float cos_theta_h)
{
  if(p == 0) // R
    return fresnel(n, cos_theta_h);

  // TT and TRT
  const float cosr = sqrtf(1.0-h*h);
  const float sin_gamma_t = h/np1;
  const float f = fresnel(n, cosr * cos_theta_d);
  const float cos2_gamma_t = 1.0f-sin_gamma_t*sin_gamma_t;
  const float cost2 = 1.0f - 1.0f/(n*n) * (1.0f - cosr*cosr);
  const float cost = sqrtf(fmaxf(0.0f, cost2));
  // [d'Eon et al. 2011] eq (13) (14) and text after that, using cos(2 gamma_t) = 2cos^2(gamma_t)-1
  const float T = (mu_a < 1e-8f) ? 1.0f : common_fasterexp(-4.0f * mu_a/cost * cos2_gamma_t);
  float A = (1-f)*(1-f) * T;
  if(p == 2) A *= f * T;
  return A;
}

#if 0
// sum from [d'Eon et al. 2011] eq (15).
// going up to iteration 40 makes it work precisely, but is still not quite as
// fast as the double bessel from the fortran library below.
static inline float bessel_I0(const float x)
{
  double sum = 1.0f;
  double x2 = x*x, xp2k = 1.0f;
  double fourpk = 1.0f;
  double fakk = 1.0f;
  for(int k=1;k<40;k++)
  {
    xp2k *= x2;
    fourpk *= 4.0f;
    fakk *= k;
    sum += xp2k / (fourpk * fakk*fakk);
  }
  return sum;
}
#endif

// From Numath Library By Tuan Dang Trong in Fortran 77.
// single-precision floating point version seems to be enough for us.
float bessel_I0(const float X)
{
  const float P1=1.0, P2=3.5156229, P3=3.0899424, P4=1.2067492;
  const float P5=0.2659732, P6=0.360768e-1, P7=0.45813e-2;
  const float Q1=0.39894228, Q2=0.1328592e-1, Q3=0.225319e-2;
  const float Q4=-0.157565e-2, Q5=0.916281e-2, Q6=-0.2057706e-1;
  const float Q7=0.2635537e-1, Q8=-0.1647633e-1, Q9=0.392377e-2;
  if (fabsf(X) < 3.75)
  {
    const float Y=(X/3.75)*(X/3.75);
    return P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))));
  }
  else
  {
    float AX=fabsf(X);
    const float Y=3.75/AX;
    const float BX=common_fasterexp(AX)/sqrtf(AX);
    AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
    return AX*BX;
  }
}

static inline float log_bessel_I0(const float x)
{
  if(x > 12.0f)
    return x + 0.5f*(-logf(2.0f*M_PI) + logf(1.0f/x) + 1.0f/(8.0f*x));
  else
    return logf(bessel_I0(x));
}


// longitudinal scattering function (theta_i, theta_o, in the plane of incoming ray and the fiber)
// used for all lobes R, TT, TRT, but with different roughness.
static inline float M(
    const float v,       // variance or roughness
    const float theta_c, // cone theta
    const float theta_o) // outgoing theta
{
  // TODO: move outside?
  float sin_theta_c, cos_theta_c;
  float sin_theta_o, cos_theta_o;
  common_sincosf(theta_c, &sin_theta_c, &cos_theta_c);
  common_sincosf(theta_o, &sin_theta_o, &cos_theta_o);
  if(v < 0.1f)
  { // special case from publons.com discussion:
    const float a = cos_theta_c * cos_theta_o / v;
    const float b = sin_theta_c * sin_theta_o / v;
    return common_fasterexp(log_bessel_I0(a) + b - 1.0f/v + 0.6931f + logf(1.0f/(2.0f*v)));
  }
  else
  {
    const float csch = 1.0f/sinhf(1.0f/v);
    const float ex = common_fasterexp(sin_theta_c*sin_theta_o/v);
    const float bessel = bessel_I0(cos_theta_c*cos_theta_o/v);
    return csch/(2.0f*v) * ex * bessel;
  }
}

// importance sample M cos^2 theta_o [d'Eon et al. 2013] eq (6) and (7)
// returns sin(theta_o)
static inline float sample_M(
    const float v,          // roughness
    const float theta_c,    // center of the cone
    const float rand1,      // random number [0,1)
    const float rand2)      // random number [0,1)
{
  float cos_th, sin_th;
  if(rand1 < 1e-4f) return 1.0f;
  common_sincosf(M_PI/2.0 - theta_c, &sin_th, &cos_th);
  // const float u = v * logf(expf(-1.0f/v) + 2.0f * rand1 * sinhf(1.0f/v));
  // numerically stable variant of the above [Jakob 2012]
  const float u = 1.0f + v*(logf(rand1) + logf(1.0f - (rand1-1.0f)/rand1 * common_fasterexp(-2.0f/v)));
  // return sin(theta_o)
  return u * cos_th + sqrtf(fmaxf(0.0f, 1.0-u*u)) * cosf(2.0f*M_PI*rand2) * sin_th;
  // const float cos2_theta_o = 1.0f - *sin_theta_o**sin_theta_o;
  // weight is 1.0f/cos2_theta_o;
}


// helper function, main location of the outgoing lobe by bravais analysis
static inline float Phi(
    const int p,     // number of lobe: 0 - R, 1 - TT, 2 - TRT
    const float h,   // offset from center of fiber, -1..1
    const float np1) // bravais equivalent to ior, perpendicular component
{
  const float gamma_i = asinf(h);
  const float gamma_t = asinf(CLAMP(h/np1, -1.0, 1.0));
  const float offset = (p == 1) ? M_PI : 0;
  return 2.0f * p * gamma_t - 2.0f * gamma_i + offset;
}

// gaussian detector [d'Eon et al 2011]
static inline float D(
    const float v,   // variance = beta^2
    const float phi)
{
  // sum from k=-inf to +inf
  // sum from k=1..inf and k=1..-inf and stop when exp is zero
  const float norm = 1.0f/sqrtf(2.0f*M_PI*v);
  float sum = 0.0f;
  const int maxterms = 10;
  for(int k=0;k<maxterms;k++)
  {
    const float t = phi + 2.0f*M_PI*k;
    const float add = common_fasterexp(-t*t/(2.0f*v));
    sum += add;
    if(add < 1e-8f) break;
  }
  for(int k=1;k<maxterms;k++)
  {
    const float t = phi - 2.0f*M_PI*k;
    const float add = common_fasterexp(-t*t/(2.0f*v));
    sum += add;
    if(add < 1e-8f) break;
  }
  return sum*norm;
}

static inline float N_p_spec(
    const int p,           // lobe: 0 - R, 1 - TT, 2 - TRT
    const float v,         // variance of gaussian or roughness of lobe
    const float phi,
    const float n,         // ior
    const float np1,       // eta', modified index of refraction
    const float np1_s,     // eta', modified index of refraction, used to evaluate A() (different for pdf computation)
    const float mu_a,      // internal absorption
    const float cos_theta_d,
    const float cos_theta_h)
{
  if(p == 0)
  { // R, closed form
    return 0.25f * fabsf(cosf(phi*.5f))/(2.0f*M_PI*M_PI);
  }
#if 0 // buys a lot of speed, but seems to be broken (some modulo madness, images look off)
  else if(p == 1)
  { // TT, closed form, more complicated
    const float a = 1.0/np1;
    // h_p looks good, but something else stays fishy here.
    const float h_p = cosf(phi*.5f)/sqrtf(1.0 + a*a - 2.0f*a*fabsf(sinf(phi*.5f)));
    const float sum_A = 
      A(0, h_p, n, np1_s, mu_a, cos_theta_d, cos_theta_h)+
      A(1, h_p, n, np1_s, mu_a, cos_theta_d, cos_theta_h)+
      A(2, h_p, n, np1_s, mu_a, cos_theta_d, cos_theta_h);
    return .5f * A(1, h_p, n, np1_s, mu_a, cos_theta_d, cos_theta_h) / sum_A;
  }
#endif
  // TRT, fall back to quadrature :(
  else
  {
    const int num = 35;
    const float dh = 1.0f/(num-1.0f);
    float sum = 0.0f;
    for(int i=0;i<num;i++)
    {
      const float h = i/(num-1.0f);
      const float Phi_ = Phi(p, h, np1);
      const float phi2 = fabsf(phi - Phi_);
      const float phi3 = fabsf(phi + Phi_);
      const float sum_A = 
        A(0, h, n, np1_s, mu_a, cos_theta_d, cos_theta_h) +
        A(1, h, n, np1_s, mu_a, cos_theta_d, cos_theta_h) +
        A(2, h, n, np1_s, mu_a, cos_theta_d, cos_theta_h);
      sum += A(p, h, n, np1_s, mu_a, cos_theta_d, cos_theta_h) * (D(v, phi2) + D(v, phi3)) / sum_A;
    }
    return .5f * sum * dh;
  }
}


// azimuthal scattering function (phi difference, orthogonal to the fiber) [d'Eon et al. 2011]
static inline float N_p(
    const int p,           // lobe: 0 - R, 1 - TT, 2 - TRT
    const float v,         // variance of gaussian or roughness of lobe
    const float phi,       // azimuthal difference angle
    const float n,         // ior
    const float np1,       // eta', modified index of refraction
    const float np1_s,     // eta', modified index of refraction, used to evaluate A() (different for pdf computation)
    const float mu_a,      // internal absorption
    const float cos_theta_d,
    const float cos_theta_h)
{
  // quadrature of
  // 0.5f * integral h=-1..1 A(p, h) * D(v_p, phi - Phi(p, h, np1)) dh
  // integrate from 0..1 instead and use symmetry
  const int num = 35;
  const float dh = 1.0f/(num-1.0f);
  float sum = 0.0f;
  for(int i=0;i<num;i++)
  {
    // const float h = -1.0f + 2.0f * i/(num-1.0f);
    const float h = i/(num-1.0f);
    const float Phi_ = Phi(p, h, np1);
    const float phi2 = fabsf(phi - Phi_);
    const float phi3 = fabsf(phi + Phi_);
    sum += A(p, h, n, np1_s, mu_a, cos_theta_d, cos_theta_h) * (D(v, phi2) + D(v, phi3));
    // sum += A(p, h, n, np1, mu_a, cos_theta_d, cos_theta_h) * D(v, phi2);
    // sum += D(v, phi2);
  }
  return .5f * sum * dh;
  // return sum * dh;
}

// sample azimuthal N function for lobe p.
// the sample weight is A(p, h) <= 1.
// returns relative phi angle
static inline float sample_N(
    const int p,              // lobe number
    const float np1,          // modified index of refraction eta'
    const float v,            // roughness
    const float h,            // pre-sampled offset from fiber centre
    const float rand_gauss)   // gaussian random number
{
  // sample phi by A() eq (8) in [d'Eon et al. 2013]
  return Phi(p, h, np1) + rand_gauss * sqrtf(v);
}

// absorption coefficients [1/mm] for melanin concentrations
// [Donner and Jensen 2006 "A Spectral BSSRDF for Shading Human Skin"]
// we need the relative absorption coefficient wrt to fiber
// width [Marschner et al. 2003], so these numbers are off by at least
// two orders of magnitude (hair: 0.1mm and we compute with 1.0f=1dm)
static inline mf_t eumelanin(mf_t lambda)
{
  float res[MF_COUNT];
  const float *lf = (float *)&lambda;
  for(int l=0;l<MF_COUNT;l++)
    // return 6.6e10f * powf(lambda, -3.33f); // 1/mm, original from paper
    res[l] = 6.6e8f * powf(lf[l], -3.33f); // adjusted relative to 0.1mm hair width
  return mf_loadu(res);
}

static inline mf_t pheomelanin(mf_t lambda)
{
  float res[MF_COUNT];
  const float *lf = (float *)&lambda;
  for(int l=0;l<MF_COUNT;l++)
    // return 2.9e14f * powf(lambda, -4.75f); // 1/mm
    res[l] = 2.9e12f * powf(lf[l], -4.75f);
  return mf_loadu(res);
}

// split a random number in two. we need moar :)
static inline void split_rand(
    const float rand,
    float *r1,
    float *r2)
{
  // odd bits
  uint32_t mantissa = touint(rand + 1.0f) & 0x7fffff;
  uint32_t x = mantissa & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff;
  *r1 = tofloat(0x3f800000 | (x<<12)) - 1.0f;

  // even bits
  x = mantissa & 0xaaaaaaaa;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff;
  *r2 = tofloat(0x3f800000 | (x<<12)) - 1.0f;
}


#ifndef HAIR_UNIT_TEST
// ===========================================================================
// public interface functions below.
// ===========================================================================

// update vertex with shading information
float prepare(path_t *p, int v, void *data)
{ 
  hair_t *hair = (hair_t *)data;
  const float n_d = 1.55f;
  const float V_d = 40.0f;
  // set volume properties on surface. we don't ever transmit, but we will use this for fresnel:
  // float r0, r1;
  // geo_line_get_radii(rt->prims, p->v[v].hit.prim, &r0, &r1);
  // const float r = p->v[v].hit.u * r1 + (1.0f-p->v[v].hit.u) * r0;
  const float t = p->v[v].hit.r;
  const float blend = 1. + .3*(.7-t*t*t*(t*(t*6.0f-15.0f)+10.0f));
  p->v[v].interior.mu_t = 
    mf_mul(mf_set1(blend), mf_fma(mf_set1(hair->pheomelanin), pheomelanin(p->lambda),
                           mf_mul(mf_set1(hair->eumelanin),   eumelanin(p->lambda))));
  float uvw[3] = {p->v[v].hit.s, p->v[v].hit.t, p->v[v].hit.r};
  const float noise = pnoise(uvw, 2, 1e-3f);
  p->v[v].interior.ior = mf_add(mf_set1(.1 * noise), spectrum_eta_from_abbe(n_d, V_d, p->lambda));
  p->v[v].shading.roughness = 0.23 * noise + 15.0f * M_PI/180.0f; // width of R lobe in degrees
  p->v[v].material_modes = s_fiber | s_glossy;
  return 1.0f;
}


// evaluate pdf of sampling
mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(p->v[v].flags & s_inside) return 0.0f; // intersection problem.
  if(p->v[v].mode != (s_fiber | s_glossy)) return 0.0f;

  float S = 0.0f;
  const float *fiber = p->v[v].hit.n;
  // take care of swapped e1, e2
  float wi[3];
  float wo[3];
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

  // get coordinate system such that dot(wi, normal) < 0 and dot(wo, normal) > 0
  float ortho[3], normal[3];
  crossproduct(wi, fiber, ortho);
  normalise(ortho);
  crossproduct(ortho, fiber, normal);
  normalise(normal);

  // incoming/outgoing theta and half and difference angles
  const float sin_theta_i = -dotproduct(wi, fiber);
  const float cos_theta_i = sqrtf(1.0f-sin_theta_i*sin_theta_i);
  const float theta_i = asinf(sin_theta_i);
  const float sin_theta_o = dotproduct(wo, fiber);
  const float theta_o = asinf(sin_theta_o);
  const float theta_d = (theta_o - theta_i)*0.5f;
  float sin_theta_d, cos_theta_d;
  common_sincosf(theta_d, &sin_theta_d, &cos_theta_d);

  // difference angle phi (computed as such with wi will be 0)
  const float sin_phi_scaled = dotproduct(ortho, wo);
  const float cos_phi_scaled = dotproduct(normal, wo);
  const float phi = atan2f(sin_phi_scaled, cos_phi_scaled);

  // ior
  const mf_t eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(mf(eta_ratio, 0) <= 0.0f) return 0.0f; // broken volume nesting
  const float n = 1./eta_ratio; 
  // bravais effective ior for refraction and perpendicular attenuation component
  const mf_t np1 = mf_sqrt(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;
  const mf_t np1_spec = mf_sqrt(n*n - sin_theta_i*sin_theta_i)/cos_theta_i;

  // roughness, tilt and internal absorption
  const float beta = p->v[v].shading.roughness;
  const float var = beta*beta;
  const float tilt = 2.0 * M_PI/180.0f;
  const mf_t mu_a = p->v[v].interior.mu_t;

  // [d'Eon et al. 2013] Sec. 3.5: evaluate with A(p,h) replaced by w_p (specular A / sum A)
  S += N_p_spec(0,     beta, phi, n, np1, np1_spec, mu_a, cos_theta_i, 1.0f) * M(var,       -theta_i+2.f*tilt, theta_o); // R
  S += N_p_spec(1, .5f*beta, phi, n, np1, np1_spec, mu_a, cos_theta_i, 1.0f) * M(var*0.25f, -theta_i-    tilt, theta_o); // TT
  S += N_p_spec(2, 2.f*beta, phi, n, np1, np1_spec, mu_a, cos_theta_i, 1.0f) * M(var*4.0f,  -theta_i-4.f*tilt, theta_o); // TRT

  // divide out jacobian from theta/phi space to solid angle
  // (this is the same cosine that will get multiplied when evaluating a G term)
  return S / sqrtf(1.0f-sin_theta_o*sin_theta_o);
}

// sample new outgoing direction, return weight
mf_t sample(path_t *path, void *data)
{
  const int v = path->length-1;
  if(path->v[v].flags & s_inside) return mf_set1(0.0f); // intersection problem.
  // we get three random numbers, but we need:
  // 1 - lobe selection
  // 1 - fiber offset
  // 2 - gaussian box-muller for 1 gaussian in N
  // 2 - M
  float rand_M1, rand_M2, rand_N1, rand_N2, rand_g1, rand_g2, rand_lobe, rand_h;
  split_rand(pointsampler(path, s_dim_scatter_mode), &rand_lobe, &rand_h);
  split_rand(pointsampler(path, s_dim_omega_x), &rand_M1, &rand_M2);
  split_rand(pointsampler(path, s_dim_omega_y), &rand_N1, &rand_N2);
  sample_gaussian(rand_N1, rand_N2, &rand_g1, &rand_g2);

  // init coordinate frame and angles:
  const float *fiber = path->v[v].hit.n;
  const float *wi = path->e[v].omega;
  float *wo = path->e[v+1].omega;

  // get coordinate system such that dot(wi, normal) < 0 and dot(wo, normal) > 0
  float ortho[3], normal[3];
  crossproduct(wi, fiber, ortho);
  normalise(ortho);
  crossproduct(ortho, fiber, normal);
  normalise(normal); // pranoia

  // incoming/outgoing theta and half and difference angles
  const float sin_theta_i = CLAMP(-dotproduct(wi, fiber), -1, 1);
  const float cos_theta_i = sqrtf(fmaxf(0.0f, 1.0f - sin_theta_i*sin_theta_i));
  const float theta_i = atan2f(sin_theta_i, cos_theta_i);

  // ior
  const mf_t eta_ratio = path_eta_ratio(path, v); // = n1/n2;
  if(mf(eta_ratio, 0) <= 0.0) return mf_set1(0.0f); // broken volume nesting
  const float n = 1.0f/mf(eta_ratio, 0);

  // roughness, tilt and internal absorption
  float beta = path->v[v].shading.roughness;
  const float tilt = 2.0 * M_PI/180.0f;
  const mf_t mu_a = path->v[v].interior.mu_t;

  // select fiber offset h
  const float h = 1.0f - 2.0f*rand_h;

  // select lobe, compute cumulative distribution function by assuming smooth scattering
  float p_cdf[3];
  const float np1_spec = sqrtf(n*n - sin_theta_i*sin_theta_i)/cos_theta_i;
  for(int k=0;k<3;k++)
    p_cdf[k] = A(k, h, n, np1_spec, mu_a, cos_theta_i, 1.0f);
  for(int k=1;k<3;k++)
    p_cdf[k] += p_cdf[k-1];
  // if(!(p_cdf[2] > 0.0f)) p_cdf[0] = p_cdf[1] = p_cdf[2] = 1.0;//return 0.0f;
  if(!(p_cdf[2] > 0.0f)) return mf_set1(0.0f);
  for(int k=0;k<2;k++)
    p_cdf[k] /= p_cdf[2];
  p_cdf[2] = 1.0f;
  const int p = sample_cdf(p_cdf, 3, rand_lobe);
  if(p == 1) beta *= 0.5f;
  else if(p == 2) beta *= 2.0f;
  const float w_p = p ? (p_cdf[p] - p_cdf[p-1]) : p_cdf[0]; // lobe selection probability

  // sample longitudinal M function:
  const float var = beta*beta;
  float theta_c = -theta_i + 2.0f*tilt;
  if(p == 1) theta_c = -theta_i-tilt;
  else if(p == 2) theta_c = -theta_i-4.0f*tilt;
  const float sin_theta_o = sample_M(var, theta_c, rand_M1, rand_M2);

  // bravais effective ior for refraction and perpendicular attenuation component
  const float theta_o = asinf(sin_theta_o);
  const float theta_d = (theta_o - theta_i)*.5f;
  float sin_theta_d, cos_theta_d;
  common_sincosf(theta_d, &sin_theta_d, &cos_theta_d);
  const float np1 = sqrtf(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;

  // sample azimuthal N function to get relative phi
  const float phi = sample_N(p, np1, beta, h, rand_g1);

  float sin_phi, cos_phi;
  common_sincosf(phi, &sin_phi, &cos_phi);
  const float cos_theta_o = sqrtf(1.0f - sin_theta_o*sin_theta_o);
  for(int k=0;k<3;k++)
    wo[k] = sin_theta_o * fiber[k] + cos_theta_o * (sin_phi * ortho[k] + cos_phi * normal[k]);

  // return clamped weight:
  const float theta_h = (theta_o + theta_i)*.5f;
  const float cos_theta_h = cosf(theta_h);
  // TODO: this needs to be spectral:
  const float weight = fminf(2.0f, A(p, h, n, np1, mu_a, cos_theta_d, cos_theta_h)/w_p);

  path->v[v].mode = s_fiber | s_glossy;

  path->v[v+1].pdf = pdf(path, v, v, v+1, data);
  return weight;
}

// evaluate curve scattering function
mf_t brdf(path_t *p, int v, void *data)
{
  if(p->v[v].flags & s_inside) return mf_set1(0.0f); // intersection problem.
  float S = 0.0f;
  const float *fiber = p->v[v].hit.n;
  const float *wi = p->e[v].omega, *wo = p->e[v+1].omega;

  // get coordinate system such that dot(wi, normal) < 0 and dot(wo, normal) > 0
  float ortho[3], normal[3];
  crossproduct(wi, fiber, ortho);
  normalise(ortho);
  crossproduct(ortho, fiber, normal);
  normalise(normal);

  // incoming/outgoing theta and half and difference angles
  const float sin_theta_i = -dotproduct(wi, fiber);
  // const float cos_theta_i = sqrtf(1.0f - sin_theta_i*sin_theta_i);
  const float theta_i = asinf(sin_theta_i);//atan2f(sin_theta_i, cos_theta_i);
  const float sin_theta_o = dotproduct(wo, fiber);
  // const float cos_theta_o = sqrtf(1.0f - sin_theta_o*sin_theta_o);
  const float theta_o = asinf(sin_theta_o);//atan2f(sin_theta_o, cos_theta_o);
  const float theta_h = (theta_o + theta_i)*0.5f;
  const float cos_theta_h = cosf(theta_h);
  const float theta_d = (theta_o - theta_i)*0.5f;
  float sin_theta_d, cos_theta_d;
  common_sincosf(theta_d, &sin_theta_d, &cos_theta_d);

  // difference angle phi (computed as such with wi will be 0)
  const float sin_phi_scaled = dotproduct(ortho, wo);
  const float cos_phi_scaled = dotproduct(normal, wo);
  const float phi = atan2f(sin_phi_scaled, cos_phi_scaled);
  // fprintf(stderr, "theta %g %g phi %g\n", theta_i, theta_o, phi);

  // ior
  const mf_t eta_ratio = path_eta_ratio(p, v); // = n1/n2;
  if(mf(eta_ratio, 0) <= 0.0) return mf_set1(0.0f); // broken volume nesting
  const mf_t n = mf_div(mf_set1(1.0f), eta_ratio);
  // bravais effective ior for refraction and perpendicular attenuation component
  const float np1 = sqrtf(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;
  // bravais effective ior for parallel attenuation component:
  //const float np2 = n*n/np1;

  // roughness, tilt and internal absorption
  const float beta = p->v[v].shading.roughness;
  const float var = beta*beta;
  const float tilt = 2.0 * M_PI/180.0f;
  const mf_t mu_a = p->v[v].interior.mu_t;

  S += N_p(0, beta,      phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var,       -theta_i+2.f*tilt, theta_o); // R
  S += N_p(1, beta*0.5f, phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var*0.25f, -theta_i-    tilt, theta_o); // TT
  S += N_p(2, beta*2.0f, phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var*4.0f,  -theta_i-4.f*tilt, theta_o); // TRT

  p->v[v].mode = s_fiber | s_glossy;
  // divide out jacobian from theta/phi space to solid angle
  // (this is the same cosine that will get multiplied when evaluating a G term)
  return S / sqrtf(1.0f-sin_theta_o*sin_theta_o);
}

// load from file
int init(FILE *f, void **data)
{
  hair_t *hair = (hair_t *)malloc(sizeof(hair_t));
  *data = hair;
  // load global parameters: melanin concentrations
  // the other two come from the global shading struct (absorption spectrum: rg, roughness: that)
  // or are left constant (tilt 2 deg, ior 1.55)
  if(fscanf(f, "%f %f", &hair->eumelanin, &hair->pheomelanin) != 2)
  {
    fprintf(stderr, "[hair] warning: could not read params (usage: hair eumelanin pheomelanin)! using defaults.\n");
    hair->pheomelanin = 0.5;
    hair->eumelanin = 0.1;
  }
  int dreggn = fscanf(f, "%*[^\n]\n");
  return dreggn == -1;
}

#else

//=========================================================================
// unit test
// clang -std=c11 -march=native -DHAIR_UNIT_TEST -D_GNU_SOURCE -I ../../ -I ../../include -I ../../build  hair.c -o hair -lm
//=========================================================================

int main(int argc, char *argv[])
{
  const int wd = 128, ht = 256;
  float *buf = (float *)malloc(sizeof(float)*3*wd*ht);

  const float n = 1.55f;
  const float beta = 15.0f * M_PI/180.0f; // 15 degrees
  const float var = beta*beta;
  const float tilt = 0.0f;
  const float pheo = 0.0f, eu = 0.0f;
  const float lambda = 550.f;
  const float mu_a = pheo * pheomelanin(lambda) + eu * eumelanin(lambda);

  const float theta_i = -M_PI/4.0f;

  double start = common_time_wallclock();
  // =========================================================================
  // evaluation:
  for(int j=0;j<ht;j++)
  {
    for(int i=0;i<wd;i++)
    {
      // outgoing theta [-pi/2..pi/2]
      const float theta_o = ((j/(ht-1.0f)) - .5f)*2.0f  * M_PI/2.0f;
      // difference phi [-pi..pi]
      const float phi = ((i/(wd-1.0f)) - .5f)*2.0f  * M_PI;
      const float theta_h = (theta_o + theta_i)*0.5f;
      const float cos_theta_h = cosf(theta_h);
      const float theta_d = (theta_o - theta_i)*0.5f;
      const float cos_theta_d = cosf(theta_d);
      const float sin_theta_d = sinf(theta_d);
      // XXX typo in [d'Eon 2013], sin^2 missing!
      const float np1 = sqrtf(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;
      // cos theta_o is needed to convert bsdf eval to sampling (which is sampling f_r * cos_o)
      // jacobian from theta/phi to solid angle is cosf(theta_o), too, but
      // we're not dividing out another factor of cos(theta_o) in M(), but only do
      // it in bsdf(), which is not called here.
      const float det = cosf(theta_o);
      buf[3*(j*wd + i) + 0] = N_p(0, beta,      phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var,       -theta_i + 2.f * tilt, theta_o)*det; // R
      buf[3*(j*wd + i) + 1] = N_p(1, beta*0.5f, phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var*0.25f, -theta_i -       tilt, theta_o)*det; // TT
      buf[3*(j*wd + i) + 2] = N_p(2, beta*2.0f, phi, n, np1, np1, mu_a, cos_theta_d, cos_theta_h) * M(var*4.0f,  -theta_i - 4.f * tilt, theta_o)*det; // TRT
    }
  }
  double end = common_time_wallclock();
  fprintf(stderr, "eval in %g sec (%g per eval)\n", end-start, (end-start)/(wd*ht));

  // compute average to check energy conservation:
  double sum[3] = {0};
  double norm = 2.0f*M_PI*M_PI/(double)(wd*ht);
  for(int k=0;k<wd*ht;k++)
    for(int c=0;c<3;c++) sum[c] += (double)buf[3*k+c];
  fprintf(stderr, "energy output: %g %g %g sum %g\n", sum[0]*norm, sum[1]*norm, sum[2]*norm, (sum[0]+sum[1]+sum[2])*norm);

  FILE *f = fopen("hair-eval.pfm", "wb");
  fprintf(f, "PF\n%d %d\n-1.0\n", wd, ht);
  fwrite(buf, 3*sizeof(float), wd*ht, f);
  fclose(f);


  start = common_time_wallclock();
  // =========================================================================
  // pdf evaluation:
  for(int j=0;j<ht;j++)
  {
    for(int i=0;i<wd;i++)
    {
      // outgoing theta [-pi/2..pi/2]
      const float theta_o = ((j/(ht-1.0f)) - .5f)*2.0f  * M_PI/2.0f;
      // difference phi [-pi..pi]
      const float phi = ((i/(wd-1.0f)) - .5f)*2.0f  * M_PI;
      const float theta_h = (theta_o + theta_i)*0.5f;
      const float cos_theta_h = cosf(theta_h);
      const float theta_d = (theta_o - theta_i)*0.5f;
      const float cos_theta_d = cosf(theta_d);
      const float sin_theta_d = sinf(theta_d);
      // XXX typo in [d'Eon 2013], sin^2 missing!
      const float np1 = sqrtf(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;
      const float np1_spec = sqrtf(n*n - sinf(theta_i)*sinf(theta_i))/cosf(theta_i);
      // cos theta_o is needed to convert bsdf eval to sampling (which is sampling f_r * cos_o)
      // jacobian from theta/phi to solid angle is cosf(theta_o), too, but
      // we're not dividing out another factor of cos(theta_o) in M(), but only do
      // it in pdf(), which is not called here.
      const float det = cosf(theta_o);
      buf[3*(j*wd + i) + 0] = N_p_spec(0, beta,      phi, n, np1, np1_spec, mu_a, cos_theta_d, cos_theta_h) * M(var,       -theta_i + 2.f * tilt, theta_o)*det; // R
      buf[3*(j*wd + i) + 1] = N_p_spec(1, beta*0.5f, phi, n, np1, np1_spec, mu_a, cos_theta_d, cos_theta_h) * M(var*0.25f, -theta_i -       tilt, theta_o)*det; // TT
      buf[3*(j*wd + i) + 2] = N_p_spec(2, beta*2.0f, phi, n, np1, np1_spec, mu_a, cos_theta_d, cos_theta_h) * M(var*4.0f,  -theta_i - 4.f * tilt, theta_o)*det; // TRT
    }
  }
  end = common_time_wallclock();
  fprintf(stderr, "pdf eval in %g sec (%g per eval)\n", end-start, (end-start)/(wd*ht));

  f = fopen("hair-pdf.pfm", "wb");
  fprintf(f, "PF\n%d %d\n-1.0\n", wd, ht);
  fwrite(buf, 3*sizeof(float), wd*ht, f);
  fclose(f);

  // =========================================================================
  // sampling:
  memset(buf, 0, 3*sizeof(float)*wd*ht);

  // incoming/outgoing theta and half and difference angles
  const float sin_theta_i = sinf(theta_i);
  const float cos_theta_i = cosf(theta_i);

  start = common_time_wallclock();

  const int num_samples = 1<<22;
  const float scale = (float)wd*ht/(float)num_samples/ (2.0f*M_PI * M_PI);
  double avgw = 0.0f;
  for(int s=0;s<num_samples;s++)
  {
    for(int p=0;p<3;p++)
    {
      float rand_M1, rand_M2, rand_N1, rand_N2, rand_g1, rand_g2, rand_lobe, rand_h;
      rand_M1 = drand48();
      rand_M2 = drand48();
      rand_N1 = drand48();
      rand_N2 = drand48();
      rand_lobe = drand48();
      rand_h = drand48();
      sample_gaussian(rand_N1, rand_N2, &rand_g1, &rand_g2);

      // select fiber offset h
      const float h = 1.0f - 2.0f*rand_h;

      // sample longitudinal M function:
      float beta2 = beta;
      if(p == 1) beta2 *= 0.5f;
      else if(p == 2) beta2 *= 2.0f;

      float theta_c = -theta_i + 2.0f*tilt;
      if(p == 1) theta_c = -theta_i-tilt;
      else if(p == 2) theta_c = -theta_i-4.0f*tilt;
      const float sin_theta_o = sample_M(beta2*beta2, theta_c, rand_M1, rand_M2);

      // bravais effective ior for refraction and perpendicular attenuation component
      const float theta_o = asinf(sin_theta_o);
      const float theta_d = (theta_o - theta_i)*.5f;
      float sin_theta_d, cos_theta_d;
      common_sincosf(theta_d, &sin_theta_d, &cos_theta_d);
      const float np1 = sqrtf(n*n - sin_theta_d*sin_theta_d)/cos_theta_d;

      // sample azimuthal N function to get relative phi
      const float phi_ = sample_N(p, np1, beta2, h, rand_g1); // will wrap around outside period a lot (we don't care for path space sampling)
      const float phi = fmodf(phi_ + M_PI + 1000.0* 2.0f*M_PI, 2.0f*M_PI) - M_PI;

      const int i = fminf(wd-1, fmaxf(0.0, (phi/(2.0f*M_PI) + .5f)*wd));
      const int j = fminf(ht-1, fmaxf(0.0, (theta_o/M_PI + .5f)*ht));

      // return clamped weight:
      const float theta_h = (theta_o + theta_i)*.5f;
      const float cos_theta_h = cosf(theta_h);
      const float weight = A(p, h, n, np1, mu_a, cos_theta_d, cos_theta_h);

      const float cos_theta_o = cosf(theta_o);
      buf[3*(wd*j + i) + p] += scale * fminf(2.0f, weight);
      avgw += fminf(2.0f, weight);
    }
  }
  end = common_time_wallclock();
  fprintf(stderr, "sample in %g sec (%g per sample)\n", end-start, (end-start)/num_samples);
  fprintf(stderr, "average weight: %g\n", avgw/num_samples);

  f = fopen("hair-sample.pfm", "wb");
  fprintf(f, "PF\n%d %d\n-1.0\n", wd, ht);
  fwrite(buf, 3*sizeof(float), wd*ht, f);
  fclose(f);


  free(buf);
  exit(0);
}
#endif

