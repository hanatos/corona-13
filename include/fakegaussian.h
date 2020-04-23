#pragma once
// fake gaussian distribution, normalised to (-inf,inf) or, equivalently,
// [-4,4]. it also samples this domain. close to a real gaussian with
// sigma=0.833 (i.e. sqrt(log(2))).

static inline float fakegaussian_pdf(const float x)
{
  const float norm_c1 = 2.0f * 0x1.1903a6p+0;
  const int i1 = 0x3f800000u, i2 = 0x40000000u;
  const int k0 = i1 - x*x * (i2 - i1);
  const int k = k0 > 0 ? k0 : 0;
  return (*(const float *)&k)/norm_c1;
}

static inline float fakegaussian_sample(float xi)
{
  float sign = 1.0f;
  if(xi >= 0.5f)
  {
    sign = -1.0f;
    xi = 2.0f*(xi-0.5f);
  }
  else xi *= 2.0f;

  const int i1 = 0x3f800000u, i2 = 0x40000000u;
  const float norm_c1 = 0x1.1903a6p+0;
  static const float cdf_lut[17] = {
    0x0p+0, 0x1.84aff2p-1, 0x1.ce84acp-1, 0x1.eaa068p-1,
    0x1.f66fcep-1, 0x1.fba148p-1, 0x1.fdf994p-1, 0x1.ff0d6p-1,
    0x1.ff8da8p-1, 0x1.ffc9ep-1, 0x1.ffe658p-1, 0x1.fff3ep-1,
    0x1.fffa56p-1, 0x1.fffd7p-1, 0x1.fffeeep-1, 0x1.ffffa6p-1,
    0x1p+0,
  };
  static const float sqrti_lut[17] = {
    0x0p+0, 0x1p+0, 0x1.6a09e6p+0, 0x1.bb67aep+0, 0x1p+1,
    0x1.1e377ap+1, 0x1.3988e2p+1, 0x1.52a7fap+1, 0x1.6a09e6p+1,
    0x1.8p+1, 0x1.94c584p+1, 0x1.a8872ap+1, 0x1.bb67aep+1,
    0x1.cd82b4p+1, 0x1.deeea2p+1, 0x1.efbdecp+1, 0x1p+2,
  };
  static const float slope_lut[17] = {
    0x1.513794p+0, 0x1.6fad2ap+1, 0x1.7286d6p+2, 0x1.73b866p+3,
    0x1.74614p+4, 0x1.74cc88p+5, 0x1.7516c2p+6, 0x1.754d32p+7,
    0x1.7576d2p+8, 0x1.7597bp+9, 0x1.75b24ep+10, 0x1.75c84ap+11,
    0x1.75dac4p+12, 0x1.75ea8p+13, 0x1.75f814p+14, 0x1.7603e6p+15,
    0x1.760e48p+16 };
  static const float a3_lut[17] = {
    -0x1.555556p-3, -0x1.555556p-4, -0x1.555556p-5, -0x1.555556p-6,
    -0x1.555556p-7, -0x1.555556p-8, -0x1.555556p-9, -0x1.555556p-10,
    -0x1.555556p-11, -0x1.555556p-12, -0x1.555556p-13, -0x1.555556p-14,
    -0x1.555556p-15, -0x1.555556p-16, -0x1.555556p-17, -0x1.555556p-18,
    -0x1.555556p-19
  };
  static const float b_lut[17] = {
    0x1p+0, 0x1.8p-1, 0x1p-1, 0x1.4p-2,
    0x1.8p-3, 0x1.cp-4, 0x1p-4, 0x1.2p-5,
    0x1.4p-6, 0x1.6p-7, 0x1.8p-8, 0x1.ap-9,
    0x1.cp-10, 0x1.ep-11, 0x1p-11, 0x1.1p-12,
    0x1.2p-13,
  };

  int i=0;
  for(;i<16;i++) if(__builtin_expect(cdf_lut[i+1] >= xi, 1)) break;

  const float cdfm = cdf_lut[i];
  const float sqrti  = sqrti_lut[i];
  const float x0 = sqrti + (xi-cdfm) * slope_lut[i];
  if(i > 3) return sign*x0; // sometimes shaves off a couple % run time, but not much.

  const float a = a3_lut[i];
  const float b = b_lut[i];
  const int k0 = i1 - x0*x0 * (i2 - i1);
  const float pdf_x0 = 1.0f/ *(const float *)&k0;
  const float bound = (cdfm - xi)*norm_c1 - (a*i + b)*sqrti;
  const float x1 = x0
    - (a*x0*x0*x0 + b*x0 + bound)*pdf_x0;
  if(i) return sign*x1;
  const float x2 = x1
    - (a*x1*x1*x1 + b*x1 + bound)*pdf_x0;
  return sign*x2;
}

