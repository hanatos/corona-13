#pragma once

#include "corona_common.h"
#include "rgb2spec.h"
#include "mf.h"

#include <math.h>

typedef enum vol_shader_id_t
{
  s_vol_shader_blackbody = 0,
  s_vol_shader_biolum    = 1,
  s_vol_shader_const     = 2,
}
vol_shader_id_t;

// shader to return radiant intensity per wavelength
// TODO: two more versions: one to return avg over 16 wavelengths (hierarchy creation)
// TODO: the other to return 8 voxels directly (for sampling)
typedef mf_t (*vol_emission_shader_t)(const float density, const float temperature, const mf_t lambda);

// TODO: identify shader by something (name, hash, 16-bin samples)

mf_t vol_shader_blackbody(const float unused_density, const float temperature, const mf_t lambda)
{
  // const double h = 6.62606957e-34; // Planck's constant [J s]
  // const double c = 299792458.0;    // speed of light [m/s]
  // const double k = 1.3807e-23;     // Boltzmann's constant [J/K]
  // const double lambda_m = lambda*1e-9; // lambda [m]
  // const double lambda2 = lambda_m*lambda_m;
  // const double lambda5 = lambda2*lambda_m*lambda2;
  // const double c1 = 2. * h * c * c / lambda5;
  // const double c2 = h * c / (lambda_m * temperature * k);
  // // convert to spectral radiance in [W/m^2 / sr / nm]
  // return c1 / (exp(c2)-1.0) * 1e-9;
  // //return c1 / (common_fasterexp(c2)-1.0) * 1e-9;
  const double h = 6.62606957e-34; // Planck's constant [J s]
  const double c = 299792458.0;    // speed of light [m/s]
  const double k = 1.3807e-23;     // Boltzmann's constant [J/K]
  const mf_t lambda2 = mf_mul(lambda, lambda);
  const mf_t lambda5 = mf_mul(lambda2, mf_mul(lambda, lambda2));
  const mf_t c1 = mf_div(mf_set1(1e45 * h * c * c), lambda5); // 1e45 to convert to nm
  const mf_t c2 = mf_div(mf_set1(h * c * 1e9 / k), mf_mul(lambda, mf_set1(temperature)));
  // convert to spectral radiance in [W/m^2 / sr / nm]
  return mf_div(c1, mf_sub(mf_exp(c2), mf_set1(1.0f))) * 1e-9f;
}

mf_t vol_shader_biolum(const float density, const float temperature, const mf_t lambda)
{
  static const float coeff[] = { 0.0, 0.5, 0.8 };
  return mf_mul(mf_set1(1000.f * temperature), mf_rgb2spec(coeff, lambda));
}

mf_t vol_shader_const(const float density, const float temperature, const mf_t lambda)
{
  return mf_set1(1.0f);
}

vol_emission_shader_t vol_shader_get(vol_shader_id_t id)
{
  switch(id)
  {
    case s_vol_shader_blackbody:
      return vol_shader_blackbody;
    case s_vol_shader_biolum:
      return vol_shader_biolum;
    case s_vol_shader_const:
      return vol_shader_const;
    default:
      return vol_shader_blackbody;
  }
}
