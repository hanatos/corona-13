#pragma once
// simple float<->half conversion routines,
// not handling special cases correctly and doing wrong rounding. whatever.

typedef union chalf_t
{
  uint16_t u;
  struct
  {
    uint16_t m : 10;
    uint16_t e : 5;
    uint16_t s : 1;
  };
}
chalf_t;

typedef union cfloat_t
{
  uint32_t u;
  float f;
  struct
  {
    uint32_t m : 23;
    uint32_t e : 8;
    uint32_t s : 1;
  };
}
cfloat_t;

// https://gist.github.com/rygorous/2156668
static inline uint16_t float_to_half(const float input)
{
  cfloat_t f32infty = { 255 << 23 };
  cfloat_t f16max = { (127 + 16) << 23 };
  cfloat_t magic = { 15 << 23 };
  cfloat_t expinf = { (255 ^ 31) << 23 };
  uint32_t sign_mask = 0x80000000u;
  chalf_t o = { 0 };
  cfloat_t f;
  f.f = input;

  uint32_t sign = f.u & sign_mask;
  f.u ^= sign;

  if (!(f.f < f32infty.u)) // Inf or NaN
    o.u = f.u ^ expinf.u;
  else
  {
    if (f.f > f16max.f) f.f = f16max.f;
    f.f *= magic.f;
  }

  o.u = f.u >> 13; // Take the mantissa bits
  o.u |= sign >> 16;
  return o.u;
}

static inline float half_to_float(const uint16_t input)
{
  static const cfloat_t magic = { 113 << 23 };
  static const uint32_t shifted_exp = 0x7c00 << 13; // exponent mask after shift
  cfloat_t o;
  chalf_t h;
  h.u = input;

  o.u = (h.u & 0x7fff) << 13;     // exponent/mantissa bits
  uint32_t exp = shifted_exp & o.u;   // just the exponent
  o.u += (127 - 15) << 23;        // exponent adjust

  // handle exponent special cases
  if (exp == shifted_exp) // Inf/NaN?
    o.u += (128 - 16) << 23;    // extra exp adjust
  else if (exp == 0) // Zero/Denormal?
  {
    o.u += 1 << 23;             // extra exp adjust
    o.f -= magic.f;             // renormalize
  }

  o.u |= (h.u & 0x8000) << 16;    // sign bit
  return o.f;
}

