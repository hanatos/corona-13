/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "points.h"
#include <emmintrin.h>
#include <inttypes.h>

#define MEXP 19937

/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/** Mersenne Exponent. The period of the sequence 
 *  is a multiple of 2^MEXP-1.
 * #define MEXP 19937 */
/** SFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define N (MEXP / 128 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define N32 (N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define N64 (N * 2)

#define POS1	122
#define SL1	18
#define SL2	1
#define SR1	11
#define SR2	1
#define MSK1	0xdfffffefU
#define MSK2	0xddfecb7fU
#define MSK3	0xbffaffffU
#define MSK4	0xbffffff6U
#define PARITY1	0x00000001U
#define PARITY2	0x00000000U
#define PARITY3	0x00000000U
#define PARITY4	0x13c9e684U


#define ALTI_SL1	{SL1, SL1, SL1, SL1}
#define ALTI_SR1	{SR1, SR1, SR1, SR1}
#define ALTI_MSK	{MSK1, MSK2, MSK3, MSK4}
#define ALTI_MSK64	{MSK2, MSK1, MSK4, MSK3}
#define ALTI_SL2_PERM	{1,2,3,23,5,6,7,0,9,10,11,4,13,14,15,8}
#define ALTI_SL2_PERM64	{1,2,3,4,5,6,7,31,9,10,11,12,13,14,15,0}
#define ALTI_SR2_PERM	{7,0,1,2,11,4,5,6,15,8,9,10,17,12,13,14}
#define ALTI_SR2_PERM64	{15,0,1,2,3,4,5,6,17,8,9,10,11,12,13,14}
#define IDSTR	"SFMT-19937:122-18-1-11-1:dfffffef-ddfecb7f-bffaffff-bffffff6"


/** 128-bit data structure */
typedef union w128_t
{
  __m128i si;
  uint32_t u[4];
}
w128_t;

typedef struct sfmt_state_t
{
  /** the 128-bit internal state array */
  w128_t sfmt[N];
  /** the 32bit integer pointer to the 128-bit internal state array */
  uint32_t *psfmt32;
#if !defined(BIG_ENDIAN64) || defined(ONLY64)
  /** the 64bit integer pointer to the 128-bit internal state array */
  uint64_t *psfmt64;
#endif
  /** index counter to the 32-bit internal state array */
  int idx;
  /** a flag: it is 0 if and only if the internal state is not yet
   * initialized. */
  int initialized;
  /** a parity check vector which certificate the period of 2^{MEXP} */
  uint32_t parity[4];
}
sfmt_state_t;

inline static int idxof(int i);
// inline static void rshift128(w128_t *out,  w128_t const *in, int shift);
// inline static void lshift128(w128_t *out,  w128_t const *in, int shift);
inline static void gen_rand_all(sfmt_state_t *s);
// inline static void gen_rand_array(sfmt_state_t *s, w128_t *array, int size);

/** 
 * @file SFMT.h 
 *
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) pseudorandom
 * number generator
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t  
#define PRIu64 "llu"
#define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#include <stdio.h>
#include <stdint.h>


/** generates a random number on [0,1)-real-interval (float) */
inline static float to_real2f(uint32_t v)
{
  v = 0x3f800000 | (v>>9); // faster than double version.
  return (*(float*)&v) - 1.0f;
  /* divided by 2^32 */
}

/**
 * This function generates and returns 32-bit pseudorandom number.
 * init_gen_rand or init_by_array must be called before this function.
 * @return 32-bit pseudorandom number
 */
static inline uint32_t gen_rand32(sfmt_state_t *s)
{
  uint32_t r;

  //assert(s->initialized);
  if (s->idx >= N32) {
    gen_rand_all(s);
    s->idx = 0;
  }
  r = s->psfmt32[s->idx++];
  return r;
}

static inline float genrand_real2f(sfmt_state_t *s)
{
  return to_real2f(gen_rand32(s));
}


/** 
 * @file  SFMT-sse2.h
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) for Intel SSE2
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * @note We assume LITTLE ENDIAN in this file
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

inline static __m128i mm_recursion(__m128i *a, __m128i *b, __m128i c,
    __m128i d, __m128i mask);

/**
 * This function represents the recursion formula.
 * @param a a 128-bit part of the interal state array
 * @param b a 128-bit part of the interal state array
 * @param c a 128-bit part of the interal state array
 * @param d a 128-bit part of the interal state array
 * @param mask 128-bit mask
 * @return output
 */
inline static __m128i mm_recursion(__m128i *a, __m128i *b, 
    __m128i c, __m128i d, __m128i mask) {
  __m128i v, x, y, z;

  x = _mm_load_si128(a);
  y = _mm_srli_epi32(*b, SR1);
  z = _mm_srli_si128(c, SR2);
  v = _mm_slli_epi32(d, SL1);
  z = _mm_xor_si128(z, x);
  z = _mm_xor_si128(z, v);
  x = _mm_slli_si128(x, SL2);
  y = _mm_and_si128(y, mask);
  z = _mm_xor_si128(z, x);
  z = _mm_xor_si128(z, y);
  return z;
}

/**
 * This function fills the internal state array with pseudorandom
 * integers.
 */
inline static void gen_rand_all(struct sfmt_state_t *s) {
  int i;
  __m128i r, r1, r2, mask;
  mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

  r1 = _mm_load_si128(&(s->sfmt[N - 2].si));
  r2 = _mm_load_si128(&(s->sfmt[N - 1].si));
  for (i = 0; i < N - POS1; i++) {
    r = mm_recursion(&(s->sfmt[i].si), &(s->sfmt[i + POS1].si), r1, r2, mask);
    _mm_store_si128(&(s->sfmt[i].si), r);
    r1 = r2;
    r2 = r;
  }
  for (; i < N; i++) {
    r = mm_recursion(&(s->sfmt[i].si), &(s->sfmt[i + POS1 - N].si), r1, r2, mask);
    _mm_store_si128(&(s->sfmt[i].si), r);
    r1 = r2;
    r2 = r;
  }
}

#if 0
/**
 * This function fills the user-specified array with pseudorandom
 * integers.
 *
 * @param array an 128-bit array to be filled by pseudorandom numbers.  
 * @param size number of 128-bit pesudorandom numbers to be generated.
 */
inline static void gen_rand_array(struct sfmt_state_t *s, w128_t *array, int size) {
  int i, j;
  __m128i r, r1, r2, mask;
  mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

  r1 = _mm_load_si128(&(s->sfmt[N - 2].si));
  r2 = _mm_load_si128(&(s->sfmt[N - 1].si));
  for (i = 0; i < N - POS1; i++) {
    r = mm_recursion(&(s->sfmt[i].si), &(s->sfmt[i + POS1].si), r1, r2, mask);
    _mm_store_si128(&array[i].si, r);
    r1 = r2;
    r2 = r;
  }
  for (; i < N; i++) {
    r = mm_recursion(&(s->sfmt[i].si), &array[i + POS1 - N].si, r1, r2, mask);
    _mm_store_si128(&array[i].si, r);
    r1 = r2;
    r2 = r;
  }
  /* main loop */
  for (; i < size - N; i++) {
    r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2,
        mask);
    _mm_store_si128(&array[i].si, r);
    r1 = r2;
    r2 = r;
  }
  for (j = 0; j < 2 * N - size; j++) {
    r = _mm_load_si128(&array[j + size - N].si);
    _mm_store_si128(&(s->sfmt[j].si), r);
  }
  for (; i < size; i++) {
    r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2,
        mask);
    _mm_store_si128(&array[i].si, r);
    _mm_store_si128(&(s->sfmt[j++].si), r);
    r1 = r2;
    r2 = r;
  }
}
#endif

/** 
 * @file  SFMT.c
 * @brief SIMD oriented Fast Mersenne Twister(SFMT)
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */
#include <string.h>
#include <assert.h>


typedef struct points_t
{
  sfmt_state_t **s;
  unsigned int num;
}
points_t;

inline static int idxof(int i) {
  return i;
}
#if 0
inline static void rshift128(w128_t *out, w128_t const *in, int shift) {
  uint64_t th, tl, oh, ol;

  th = ((uint64_t)in->u[3] << 32) | ((uint64_t)in->u[2]);
  tl = ((uint64_t)in->u[1] << 32) | ((uint64_t)in->u[0]);

  oh = th >> (shift * 8);
  ol = tl >> (shift * 8);
  ol |= th << (64 - shift * 8);
  out->u[1] = (uint32_t)(ol >> 32);
  out->u[0] = (uint32_t)ol;
  out->u[3] = (uint32_t)(oh >> 32);
  out->u[2] = (uint32_t)oh;
}
/**
 * This function simulates SIMD 128-bit left shift by the standard C.
 * The 128-bit integer given in in is shifted by (shift * 8) bits.
 * This function simulates the LITTLE ENDIAN SIMD.
 * @param out the output of this function
 * @param in the 128-bit data to be shifted
 * @param shift the shift value
 */
inline static void lshift128(w128_t *out, w128_t const *in, int shift) {
  uint64_t th, tl, oh, ol;

  th = ((uint64_t)in->u[3] << 32) | ((uint64_t)in->u[2]);
  tl = ((uint64_t)in->u[1] << 32) | ((uint64_t)in->u[0]);

  oh = th << (shift * 8);
  ol = tl << (shift * 8);
  oh |= tl >> (64 - shift * 8);
  out->u[1] = (uint32_t)(ol >> 32);
  out->u[0] = (uint32_t)ol;
  out->u[3] = (uint32_t)(oh >> 32);
  out->u[2] = (uint32_t)oh;
}
#endif

static inline void period_certification(sfmt_state_t *s)
{
  int inner = 0;
  int i, j;
  uint32_t work;

  for (i = 0; i < 4; i++)
    inner ^= s->psfmt32[idxof(i)] & s->parity[i];
  for (i = 16; i > 0; i >>= 1)
    inner ^= inner >> i;
  inner &= 1;
  /* check OK */
  if (inner == 1) {
    return;
  }
  /* check NG, and modification */
  for (i = 0; i < 4; i++) {
    work = 1;
    for (j = 0; j < 32; j++) {
      if ((work & s->parity[i]) != 0) {
        s->psfmt32[idxof(i)] ^= work;
        return;
      }
      work = work << 1;
    }
  }
}


/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 *
 * @param seed a 32-bit integer used as the seed.
 */
static inline void init_gen_rand(sfmt_state_t *s, uint32_t seed)
{
  int i;

  s->psfmt32[idxof(0)] = seed;
  for (i = 1; i < N32; i++) {
    s->psfmt32[idxof(i)] = 1812433253UL * (s->psfmt32[idxof(i - 1)] 
        ^ (s->psfmt32[idxof(i - 1)] >> 30))
      + i;
  }
  s->idx = N32;
  period_certification(s);
  s->initialized = 1;
}

/*----------------
  PUBLIC FUNCTIONS
  ----------------*/

points_t *points_init(const unsigned int num_threads, uint64_t frame)
{
  points_t *p = (points_t *)malloc(sizeof(points_t));
  // printf("[points_init] initing mersenne twister for %d threads...\n", num_threads);
  sfmt_state_t *states = (sfmt_state_t *) common_alloc(16, sizeof(sfmt_state_t)*num_threads);
  p->s = (sfmt_state_t **) malloc(sizeof(sfmt_state_t *)*num_threads);
  p->num = num_threads;

  uint64_t seed = frame;
  for(int i=0;i<(int)num_threads;i++)
  {
    p->s[i] = states + i;
    p->s[i]->psfmt32 = &(p->s[i]->sfmt[0].u[0]);
    p->s[i]->initialized = 0;
    p->s[i]->parity[0] = PARITY1;
    p->s[i]->parity[1] = PARITY2;
    p->s[i]->parity[2] = PARITY3;
    p->s[i]->parity[3] = PARITY4;
    init_gen_rand(p->s[i], seed);
    seed ^= (seed << 1)^1337;
  }
  return p;
}

void points_cleanup(points_t *p)
{
  free(p->s[0]);
  free(p->s);
  free(p);
}

float points_rand(points_t *p, const unsigned int thread_num)
{
  return genrand_real2f(p->s[thread_num]);
}

void points_print_info(FILE *fd)
{
  fprintf(fd, "points   : simple and fast mersenne twister\n");
}

void points_set_state(points_t *p, int thread_num, uint64_t s0, uint64_t s1) {}

// clean up after all this bs:
#undef MEXP
#undef N
#undef N32
#undef N64
#undef POS1
#undef SL1
#undef SL2
#undef SR1
#undef SR2
#undef MSK1
#undef MSK2
#undef MSK3
#undef MSK4	
#undef PARITY1
#undef PARITY2
#undef PARITY3
#undef PARITY4
#undef ALTI_SL1
#undef ALTI_SR1
#undef ALTI_MSK
#undef ALTI_MSK64
#undef ALTI_SL2_PERM
#undef ALTI_SL2_PERM64
#undef ALTI_SR2_PERM
#undef ALTI_SR2_PERM64
#undef IDSTR
#undef BIG_ENDIAN64
