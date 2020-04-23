#pragma once
#include <immintrin.h>

#if !defined(MF_COUNT)
#define MF_COUNT 1 // conservative default because it doesn't work completely yet
// #if defined(__AVX__)
// #define MF_COUNT 8
// #elif defined(__SSE_4_2__)
// #define MF_COUNT 4
// #else
// #define MF_COUNT 1
// #endif
#endif

// multi float value, used for hero wavelength sampling
// wrapped here so we can use scalars, avx, or sse.

// implemented using eight values and avx
// (but wrapped around so we could implement using SSE4.2 or scalar values)

// avx version, 8 wavelengths
#if MF_COUNT==8
#define mf_rgb2spec rgb2spec_eval_avx
typedef __m256 mf_t;
#define mf_rcp _mm256_rcp_ps
#define mf_set1 _mm256_set1_ps
#define mf_add _mm256_add_ps
#define mf_sub _mm256_sub_ps
#define mf_fma _mm256_fmadd_ps
#define mf_mul _mm256_mul_ps
#define mf_div _mm256_div_ps
#define mf_loadu _mm256_loadu_ps
#define mf_load _mm256_load_ps
#define mf_store _mm256_store_ps
#define mf_rsqrt _mm256_rsqrt_ps
#define mf_sqrt  _mm256_sqrt_ps
#define mf_exp   exp256_ps
#define mf_clamp(A,B,C) _mm256_min_ps(_mm256_max_ps(A, _mm256_set1_ps(B)), _mm256_set1_ps(C))
#define mf_abs(A) _mm256_and_ps(_mm256_set1_epi32(0x7fffffffu), A)
#define mf_neg(A) _mm256_xor_ps(_mm256_set1_epi32(0x80000000u), A)

#define mf_hero _mm256_set_epi32(0u,~0u,~0u,~0u,~0u,~0u,~0u,~0u)

static inline float mf_hsum(mf_t a)
{
  __m256 t1 = _mm256_hadd_ps(a,a);
  __m256 t2 = _mm256_hadd_ps(t1,t1);
  __m128 t3 = _mm256_extractf128_ps(t2,1);
  __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2),t3);
  return _mm_cvtss_f32(t4);
}

#if 0
_PS256_CONST(1  , 1.0f);
_PS256_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS256_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS256_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS256_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS256_CONST_TYPE(sign_mask, int, 0x80000000);
_PS256_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST256(0, 0);
_PI32_CONST256(1, 1);
_PI32_CONST256(inv1, ~1);
_PI32_CONST256(2, 2);
_PI32_CONST256(4, 4);
_PI32_CONST256(0x7f, 0x7f);

_PS256_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS256_CONST(cephes_log_p0, 7.0376836292E-2);
_PS256_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS256_CONST(cephes_log_p2, 1.1676998740E-1);
_PS256_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS256_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS256_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS256_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS256_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS256_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS256_CONST(cephes_log_q1, -2.12194440e-4);
_PS256_CONST(cephes_log_q2, 0.693359375);
static inline __m256 log256_ps(__m256 x) {
  __m256i imm0;
  __m256  one = _mm256_set1_ps(1.0f);

  //v8sf invalid_mask = _mm256_cmple_ps(x, _mm256_setzero_ps());
  v8sf invalid_mask = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LE_OS);

  x = _mm256_max_ps(x, *(v8sf*)_ps256_min_norm_pos);  /* cut off denormalized stuff */

  // can be done with AVX2
  imm0 = _mm256_srli_epi32(_mm256_castps_si256(x), 23);

  /* keep only the fractional part */
  x = _mm256_and_ps(x, *(v8sf*)_ps256_inv_mant_mask);
  x = _mm256_or_ps(x, *(v8sf*)_ps256_0p5);

  // this is again another AVX2 instruction
  imm0 = _mm256_sub_epi32(imm0, *(v8si*)_pi32_256_0x7f);
  v8sf e = _mm256_cvtepi32_ps(imm0);

  e = _mm256_add_ps(e, one);

  /* part2:
     if( x < SQRTHF ) {
       e -= 1;
       x = x + x - 1.0;
     } else { x = x - 1.0; }
  */
  //v8sf mask = _mm256_cmplt_ps(x, *(v8sf*)_ps256_cephes_SQRTHF);
  v8sf mask = _mm256_cmp_ps(x, *(v8sf*)_ps256_cephes_SQRTHF, _CMP_LT_OS);
  v8sf tmp = _mm256_and_ps(x, mask);
  x = _mm256_sub_ps(x, one);
  e = _mm256_sub_ps(e, _mm256_and_ps(one, mask));
  x = _mm256_add_ps(x, tmp);

  v8sf z = _mm256_mul_ps(x,x);

  v8sf y = *(v8sf*)_ps256_cephes_log_p0;
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p1);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p2);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p3);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p4);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p5);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p6);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p7);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p8);
  y = _mm256_mul_ps(y, x);

  y = _mm256_mul_ps(y, z);

  tmp = _mm256_mul_ps(e, *(v8sf*)_ps256_cephes_log_q1);
  y = _mm256_add_ps(y, tmp);


  tmp = _mm256_mul_ps(z, *(v8sf*)_ps256_0p5);
  y = _mm256_sub_ps(y, tmp);

  tmp = _mm256_mul_ps(e, *(v8sf*)_ps256_cephes_log_q2);
  x = _mm256_add_ps(x, y);
  x = _mm256_add_ps(x, tmp);
  x = _mm256_or_ps(x, invalid_mask); // negative arg will be NAN
  return x;
}
#endif

static inline __m256 exp256_ps(__m256 x) {
/* Modified code from this source: https://github.com/reyoung/avx_mathfun
   AVX implementation of exp
   Based on "sse_mathfun.h", by Julien Pommier
   http://gruntthepeon.free.fr/ssemath/
   Copyright (C) 2012 Giovanni Garberoglio
   Interdisciplinary Laboratory for Computational Science (LISC)
   Fondazione Bruno Kessler and University of Trento
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license) */

  __m256   exp_hi        = _mm256_set1_ps(88.3762626647949f);
  __m256   exp_lo        = _mm256_set1_ps(-88.3762626647949f);

  __m256   cephes_LOG2EF = _mm256_set1_ps(1.44269504088896341f);
  __m256   inv_LOG2EF    = _mm256_set1_ps(0.693147180559945f);

  __m256   cephes_exp_p0 = _mm256_set1_ps(1.9875691500E-4);
  __m256   cephes_exp_p1 = _mm256_set1_ps(1.3981999507E-3);
  __m256   cephes_exp_p2 = _mm256_set1_ps(8.3334519073E-3);
  __m256   cephes_exp_p3 = _mm256_set1_ps(4.1665795894E-2);
  __m256   cephes_exp_p4 = _mm256_set1_ps(1.6666665459E-1);
  __m256   cephes_exp_p5 = _mm256_set1_ps(5.0000001201E-1);
  __m256   fx;
  __m256i  imm0;
  __m256   one           = _mm256_set1_ps(1.0f);

  x     = _mm256_min_ps(x, exp_hi);
  x     = _mm256_max_ps(x, exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx     = _mm256_mul_ps(x, cephes_LOG2EF);
  fx     = _mm256_round_ps(fx, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);
  __m256  z      = _mm256_mul_ps(fx, inv_LOG2EF);
  x      = _mm256_sub_ps(x, z);
  z      = _mm256_mul_ps(x,x);

  __m256  y      = cephes_exp_p0;
  y      = _mm256_mul_ps(y, x);
  y      = _mm256_add_ps(y, cephes_exp_p1);
  y      = _mm256_mul_ps(y, x);
  y      = _mm256_add_ps(y, cephes_exp_p2);
  y      = _mm256_mul_ps(y, x);
  y      = _mm256_add_ps(y, cephes_exp_p3);
  y      = _mm256_mul_ps(y, x);
  y      = _mm256_add_ps(y, cephes_exp_p4);
  y      = _mm256_mul_ps(y, x);
  y      = _mm256_add_ps(y, cephes_exp_p5);
  y      = _mm256_mul_ps(y, z);
  y      = _mm256_add_ps(y, x);
  y      = _mm256_add_ps(y, one);

  /* build 2^n */
  imm0   = _mm256_cvttps_epi32(fx);
  imm0   = _mm256_add_epi32(imm0, _mm256_set1_epi32(0x7f));
  imm0   = _mm256_slli_epi32(imm0, 23);
  __m256  pow2n  = _mm256_castsi256_ps(imm0);
  y      = _mm256_mul_ps(y, pow2n);
  return y;
}

// select component
#define mf(A,B) (_mm256_extractf128_ps(A,(B)>>2))[B&3]

// comparisons:
#define mf_any(A)  (_mm256_movemask_ps(A))
#define mf_all(A)  (_mm256_movemask_ps(A)==0xff)
#define mf_gt(A,B) _mm256_cmp_ps(A, B, _CMP_GT_OS)
#define mf_lt(A,B) _mm256_cmp_ps(A, B, _CMP_LT_OS)
#define mf_eq(A,B) _mm256_cmp_ps(A, B, _CMP_EQ_OS)
#define mf_gte(A,B) _mm256_cmp_ps(A, B, _CMP_GE_OS)
#define mf_lte(A,B) _mm256_cmp_ps(A, B, _CMP_LE_OS)

// select masked based on comparison
#define mf_select(A,B,useA) _mm256_or_ps(_mm256_and_ps(useA, A),_mm256_andnot_ps(useA, B))
#define mf_or(A,B)  _mm256_or_ps(A,B)
#define mf_and(A,B) _mm256_and_ps(A,B)
#define mf_not(A)   _mm256_xor_ps(A,_mm256_set1_epi32(0xffffffffu))

// we'll need the double case, too (for pdfs and measurement contributions)
typedef struct md_t
{
  __m256d v0;
  __m256d v1;
}
md_t;

#define md_set1(A)  (md_t){.v0=_mm256_set1_pd(A),              .v1=_mm256_set1_pd(A)}
#define md_rcp(A)   (md_t){.v0=_mm256_rcp_pd((A).v0),          .v1=_mm256_rcp_pd((A).v1)}
#define md_add(A,B) (md_t){.v0=_mm256_add_pd((A).v0,(B).v0),   .v1=_mm256_add_pd((A).v1,(B).v1)}
#define md_sub(A,B) (md_t){.v0=_mm256_sub_pd((A).v0,(B).v0),   .v1=_mm256_sub_pd((A).v1,(B).v1)}
#define md_div(A,B) (md_t){.v0=_mm256_div_pd((A).v0,(B).v0),   .v1=_mm256_div_pd((A).v1,(B).v1)}
#define md_mul(A,B) (md_t){.v0=_mm256_mul_pd((A).v0,(B).v0),   .v1=_mm256_mul_pd((A).v1,(B).v1)}
#define md_sqrt(A)  (md_t){.v0=_mm256_sqrt_pd((A).v0),         .v1=_mm256_sqrt_pd((A).v1)}
#define md_rsqrt(A) (md_t){.v0=_mm256_rsqrt_pd((A).v0),        .v1=_mm256_rsqrt_pd((A).v1)}
#define md_fma(A,B,C) (md_t){.v0=_mm256_fmadd_pd((A).v0,(B).v0,(C).v0), .v1=_mm256_fmadd_pd((A).v1,(B).v1,(C).v1)}

#define md_loadu(A)  (md_t){.v0=_mm256_loadu_pd(A),   .v1=_mm256_loadu_pd((A)+8)}

#define md(A,B) ((B)<4?  _mm256_extractf128_ps((A).v0,(B)>>2)[(B)&3] : _mm256_extractf128_ps((A).v1,((B)&4)>>2)[(B)&3])

static inline double
md_hsum(md_t a)
{
  __m128d vlow  = _mm256_castpd256_pd128(a.v0);
  __m128d vhigh = _mm256_extractf128_pd(v0, 1); // high 128
  vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  double res = _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar

  vlow  = _mm256_castpd256_pd128(a.v1);
  vhigh = _mm256_extractf128_pd(v1, 1); // high 128
  vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
  high64 = _mm_unpackhi_pd(vlow, vlow);
  return res + _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}


// convert float to double
#define mf_2d(A)  (md_t){.v0=_mm256_cvtps_pd(_mm256_castps256_ps128(A)), .v1=_mm256_cvtps_pd(_mm256_extractf128_ps(A,1))}

// convert double to float
#define md_2f(A)  _mm256_insertf128_ps(_mm256_castps128_ps256( _mm256_cvtpd_ps((A).v0) ), _mm256_cvtpd_ps((A).v1), 1)

#elif MF_COUNT==4
#define mf_rgb2spec rgb2spec_eval_sse
typedef __m128 mf_t;
#define mf_set1 _mm_set1_ps
#define mf_rcp _mm_rcp_ps
#define mf_add _mm_add_ps
#define mf_sub _mm_sub_ps
#define mf_fma _mm_fmadd_ps
#define mf_mul _mm_mul_ps
#define mf_div _mm_div_ps
#define mf_load _mm_load_ps
#define mf_loadu _mm_loadu_ps
#define mf_store _mm_store_ps
#define mf_rsqrt _mm_rsqrt_ps
#define mf_sqrt _mm_sqrt_ps
#define mf_exp   exp128_ps

#define mf_clamp(A,B,C) _mm_min_ps(_mm_max_ps(A, _mm_set1_ps(B)), _mm_set1_ps(C))
#define mf_abs(A) _mm_and_ps(_mm_set1_epi32(0x7fffffffu), A)
#define mf_neg(A) _mm_xor_ps(_mm_set1_epi32(0x80000000u), A)

#define mf_hero _mm_set_epi32(0u,~0u,~0u,~0u)

static inline float mf_hsum(mf_t a)
{
  __m128 t1 = _mm_hadd_ps(a,a);
  __m128 t2 = _mm_hadd_ps(t1,t1);
  return _mm_cvtss_f32(t2);
}

static inline __m128 exp128_ps(__m128 x) {
/* Modified code from this source: https://github.com/reyoung/avx_mathfun
   AVX implementation of exp
   Based on "sse_mathfun.h", by Julien Pommier
   http://gruntthepeon.free.fr/ssemath/
   Copyright (C) 2012 Giovanni Garberoglio
   Interdisciplinary Laboratory for Computational Science (LISC)
   Fondazione Bruno Kessler and University of Trento
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license) */

  __m128   exp_hi        = _mm_set1_ps(88.3762626647949f);
  __m128   exp_lo        = _mm_set1_ps(-88.3762626647949f);

  __m128   cephes_LOG2EF = _mm_set1_ps(1.44269504088896341f);
  __m128   inv_LOG2EF    = _mm_set1_ps(0.693147180559945f);

  __m128   cephes_exp_p0 = _mm_set1_ps(1.9875691500E-4);
  __m128   cephes_exp_p1 = _mm_set1_ps(1.3981999507E-3);
  __m128   cephes_exp_p2 = _mm_set1_ps(8.3334519073E-3);
  __m128   cephes_exp_p3 = _mm_set1_ps(4.1665795894E-2);
  __m128   cephes_exp_p4 = _mm_set1_ps(1.6666665459E-1);
  __m128   cephes_exp_p5 = _mm_set1_ps(5.0000001201E-1);
  __m128   fx;
  __m128i  imm0;
  __m128   one           = _mm_set1_ps(1.0f);

  x     = _mm_min_ps(x, exp_hi);
  x     = _mm_max_ps(x, exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx     = _mm_mul_ps(x, cephes_LOG2EF);
  fx     = _mm_round_ps(fx, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);
  __m128  z      = _mm_mul_ps(fx, inv_LOG2EF);
  x      = _mm_sub_ps(x, z);
  z      = _mm_mul_ps(x,x);

  __m128  y      = cephes_exp_p0;
  y      = _mm_mul_ps(y, x);
  y      = _mm_add_ps(y, cephes_exp_p1);
  y      = _mm_mul_ps(y, x);
  y      = _mm_add_ps(y, cephes_exp_p2);
  y      = _mm_mul_ps(y, x);
  y      = _mm_add_ps(y, cephes_exp_p3);
  y      = _mm_mul_ps(y, x);
  y      = _mm_add_ps(y, cephes_exp_p4);
  y      = _mm_mul_ps(y, x);
  y      = _mm_add_ps(y, cephes_exp_p5);
  y      = _mm_mul_ps(y, z);
  y      = _mm_add_ps(y, x);
  y      = _mm_add_ps(y, one);

  /* build 2^n */
  imm0   = _mm_cvttps_epi32(fx);
  imm0   = _mm_add_epi32(imm0, _mm_set1_epi32(0x7f));
  imm0   = _mm_slli_epi32(imm0, 23);
  __m128 pow2n  = _mm_castsi128_ps(imm0);
  y      = _mm_mul_ps(y, pow2n);
  return y;
}

// select component
#define mf(A,B) ((A)[B])

// comparisons:
#define mf_any(A)  (_mm_movemask_ps(A))
#define mf_all(A)  (_mm_movemask_ps(A)==0xf)
#define mf_gt(A,B) _mm_cmpgt_ps(A, B)
#define mf_lt(A,B) _mm_cmplt_ps(A, B)
#define mf_eq(A,B) _mm_cmpeq_ps(A, B)
#define mf_gte(A,B) _mm_cmpge_ps(A, B)
#define mf_lte(A,B) _mm_cmple_ps(A, B)

// select masked based on comparison
#define mf_select(A,B,useA) _mm_or_ps(_mm_and_ps(useA, A),_mm_andnot_ps(useA, B))
#define mf_or(A,B)  _mm_or_ps(A,B)
#define mf_and(A,B) _mm_and_ps(A,B)
#define mf_not(A)   _mm_xor_ps(A,_mm_set1_epi32(0xffffffffu))

// we'll need the double case, too (for pdfs and measurement contributions)
typedef __m256d md_t;

#define md_rcp _mm256_rcp_pd
#define md_set1 _mm256_set1_pd
#define md_add _mm256_add_pd
#define md_sub _mm256_sub_pd
#define md_fma _mm256_fmadd_pd
#define md_mul _mm256_mul_pd
#define md_div _mm256_div_pd
#define md_loadu _mm256_loadu_pd
#define md_load _mm256_load_pd
#define md_store _mm256_store_pd
#define md_rsqrt _mm256_rsqrt_pd
#define md_sqrt  _mm256_sqrt_pd

#define md(A,B) _mm256_extractf128_pd(A,(B)>>1)[B&1]

static inline double md_hsum(__m256d v)
{ // https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);      // high 128
  vlow  = _mm_add_pd(vlow, vhigh);                  // reduce down to 128
  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

// convert float to double
#define mf_2d(A)  _mm256_cvtps_pd(A)

// convert double to float
#define md_2f(A)  _mm256_cvtpd_ps(A)



#elif MF_COUNT==1 // scalar
#define mf_rgb2spec rgb2spec_eval_fast
typedef float mf_t;
#define mf_rcp(A) (1.0f/(A))
#define mf_set1(A) (A)
#define mf_add(A,B) ((A)+(B))
#define mf_sub(A,B) ((A)-(B))
#define mf_fma(A,B,C)  ((A)*(B) + (C))
#define mf_mul(A,B)  ((A)*(B))
#define mf_div(A,B) ((A)/(B))
#define mf_loadu(A) (A[0])
#define mf_load(A) (A[0])
// #define mf_store ..uhm?
#define mf_rsqrt(A) (1.0f/sqrtf(A))
#define mf_sqrt(A)  sqrtf(A)
#define mf_exp(A)   expf(A)
#define mf_clamp(A,B,C) CLAMP(A, B, C)
#define mf_abs(A) fabsf(A)
#define mf_exp(A) expf(A)

#define mf_hsum(A) (A)

#define mf(A,B) (A)

// comparisons:
#define mf_any(A)   (A)
#define mf_all(A)   (A)
#define mf_gt(A,B)  ((A)> (B))
#define mf_lt(A,B)  ((A)< (B))
#define mf_eq(A,B)  ((A)==(B))
#define mf_gte(A,B) ((A)>=(B))
#define mf_lte(A,B) ((A)<=(B))
// select masked based on comparison
#define mf_select(A,B,useA) ((int)(useA) ? A : B)
#define mf_or(A,B) ((int)(A)|(int)(B))
#define mf_and(A,B) ((int)(A)&(int)(B))
#define mf_not(A) (!(int)(A))
#define mf_neg(A) (-(A))
#define mf_hero 0

typedef double md_t;
#define md_set1(A)  (A)
#define md_rcp(A)   (1.0/(A))
#define md_add(A,B) ((A)+(B))
#define md_sub(A,B) ((A)-(B))
#define md_fma(A,B) ((A)*(B)+(C))
#define md_div(A,B) ((A)/(B))
#define md_mul(A,B) ((A)*(B))

#define md_loadu(A)  (A)
#define md(A,B) (A)
#define mf_2d(A)  ((double)(A))
#define md_2f(A)  ((float)(A))

#define md_hsum(A) (A)

#else
#error "unknown MF_COUNT"
#endif
