/*
    This file is part of corona-13.
    copyright (c) 2004-2016 johannes hanika

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <fcntl.h>
#include <xmmintrin.h>

struct sampler_t;
struct view_t;
struct render_t;
struct render_tls_t;
struct accel_t;
struct points_t;
struct shader_t;
struct threads_t;
struct prims_t;
struct display_t;
struct pointsampler_t;
struct lights_t;
struct rgb2spec_t;

typedef struct primid_t
{
  uint32_t extra : 3;     // extra per-primitive bits (e.g. used to store meshing info for shells)
  uint32_t shapeid : 29;
  uint32_t vi : 28;       // points into vtxidx, which is used to find the vertex location and normals and such
  uint32_t mb : 1;        // has motion blur?
  uint32_t vcnt : 3;      // wraps the type of primitive (tri, quad, or..)
}
primid_t;

#define INVALID_PRIMID (primid_t){7u, 536870911u, 268435455u, 1u, 7u}
static inline int primid_invalid(primid_t primid)
{
#ifdef __cplusplus
  // c plus plus has a small dick:
  return primid.extra == 7 && primid.shapeid == 536870911u && primid.vi == 268435455u && primid.mb == 1u && primid.vcnt == 7u;
#else
  return (*(uint64_t *)&primid) == (*(uint64_t *)&INVALID_PRIMID);
#endif
}

static inline int primid_eq(primid_t pid1, primid_t pid2)
{
  return (*(uint64_t *)&pid1) == (*(uint64_t *)&pid2);
}

// global access to data
typedef struct rt_t
{
  int argc;
  char **argv;
  int num_threads;
  int *shader_index;
  uint64_t frames;
  uint64_t anim_frame;     // used to scramble random numbers
  uint64_t batch_frames;   // increase scaling for many cores by batching jobs
  float epsilon;

  struct shader_t *shader;
  struct prims_t *prims;
  struct accel_t *accel;
  struct sampler_t *sampler;
  struct pointsampler_t *pointsampler;
  struct view_t *view;
  struct render_t *render;
  struct points_t *points;
  struct threads_t *threads;
  struct display_t *display;
  struct lights_t *lights;
  struct rgb2spec_t *rgb2spec;

  char searchpath[256];
  char basename[256];
  char output_filename[256]; // output render file name
  char camformat[256];       // format string for multiple cameras

  // flags:
  char play;
  char quit;
  char screenshot;
}
rt_t;

// global thread local data
typedef struct rt_tls_t
{
  int tid;
  struct render_tls_t *render;
}
rt_tls_t;

extern rt_t rt;
extern __thread rt_tls_t rt_tls;

typedef struct ray_t
{
  float pos[3];
  float dir[3];
  float time;
  float min_dist;
  primid_t ignore;
}
ray_t;

typedef struct hit_t
{
  primid_t prim;    // id pointing to geometric primitive
  float u, v, w;    // texture coordinates within triangle, w is depth in shells
  float r, s, t;    // texture coordinates within texture, using uv-set of geometry (r is along hair fibers in volume)
  float a[3], b[3]; // tangent space vectors
  float n[3];       // shading normal
  float x[3];       // vertex position
  float gn[3];      // geometric normal
  int shader;       // id pointing to material
  float dist;       // ray tracing distance, meaningless for path space purposes
}
hit_t;

// minimal math:
#ifndef M_PI
  # define M_E            2.7182818284590452354   /* e */
  # define M_LOG2E        1.4426950408889634074   /* log_2 e */
  # define M_LOG10E       0.43429448190325182765  /* log_10 e */
  # define M_LN2          0.69314718055994530942  /* log_e 2 */
  # define M_LN10         2.30258509299404568402  /* log_e 10 */
  # define M_PI           3.14159265358979323846  /* pi */
  # define M_PI_2         1.57079632679489661923  /* pi/2 */
  # define M_PI_4         0.78539816339744830962  /* pi/4 */
  # define M_1_PI         0.31830988618379067154  /* 1/pi */
  # define M_2_PI         0.63661977236758134308  /* 2/pi */
  # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  # define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
  # define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

// these work for float, double, int ;)
#define crossproduct(v1, v2, res) \
  (res)[0] = (v1)[1]*(v2)[2] - (v2)[1]*(v1)[2];\
  (res)[1] = (v1)[2]*(v2)[0] - (v2)[2]*(v1)[0];\
  (res)[2] = (v1)[0]*(v2)[1] - (v2)[0]*(v1)[1]

#define dotproduct(u, v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CLAMP(a, m, M) MIN(MAX(a, m), M)

static inline void normalise(float *f)
{
  const float len = 1.0f/sqrtf(dotproduct(f, f));
  for(int k=0;k<3;k++) f[k] *= len;
}

static inline void get_perpendicular(const float* n, float* res)
{
  if(fabsf(n[1]) < 0.5)
  {
    float up[3] = {0, 1, 0};
    crossproduct(n, up, res);
  }
  else
  {
    float rg[3] = {1, 0, 0};
    crossproduct(n, rg, res);
  }
}

static inline void get_onb(const float *n, float *u, float *v)
{
  get_perpendicular(n, u);
  normalise(u);
  crossproduct(n, u, v);
  // now v should be normalised already (n and u were perpendicular and normalised)
}

static inline void get_scrambled_onb(const float scramble, const float *n, float *u, float *v)
{
  if(fabsf(n[1]) < scramble)
  {
    float up[3] = {0, 1, 0};
    crossproduct(n, up, u);
  }
  else
  {
    float rg[3] = {1, 0, 0};
    crossproduct(n, rg, u);
  }
  normalise(u);
  crossproduct(n, u, v);
  // now v should be normalised already (n and u were perpendicular and normalised)
}

// convenience struct for SSE things
typedef union
{
  __m128 m;
  float f[4];
  unsigned int i[4];
}
float4_t;

static inline float tofloat(uint32_t i)
{
  union { uint32_t i; float f; } u = {.i=i};
  return u.f;
}

static inline uint32_t touint(float f)
{
  union { float f; uint32_t i; } u = {.f=f};
  return u.i;
}

static inline double tofloat64(uint64_t i)
{
  union { uint64_t i; double f; } u = {.i=i};
  return u.f;
}

static inline uint64_t touint64(double f)
{
  union { double f; uint64_t i; } u = {.f=f};
  return u.i;
}

static inline void *common_alloc(size_t align, size_t size)
{
#ifdef __APPLE__
  return malloc(size); // mac
#else
#if 0
  void *p;
  int fail = posix_memalign(&p, align, size);
  if(fail) return 0;
  return p;
#else
  // address sanitizer requires size of aligned blocks to be aligned
  return aligned_alloc(align, (size + (align-1)) & ~(align-1));
#endif
#endif
}

static inline int common_get_threadid()
{
  return rt_tls.tid;
}

int common_load_scene(FILE *f);

void common_write_sidecar(const char *tempname);

static inline double common_time_wallclock()
{
  struct timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec - 1290608000 + (1.0/1000000.0)*time.tv_usec;
}

static inline double common_time_user()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * (1.0/1000000.0);
}

static inline void common_atomic_min(float *f, float f2)
{
  uint32_t *i = (uint32_t *)f;
  float oldf = f[0];
  uint32_t newi = touint(f2), oldi = touint(oldf);
  while(oldf > f2)
  {
    __sync_bool_compare_and_swap(i, oldi, newi);
    oldf = f[0];
    oldi = touint(oldf);
  }
}

static inline void common_atomic_max(float *f, float f2)
{
  uint32_t *i = (uint32_t *)f;
  float oldf = f[0];
  uint32_t newi = touint(f2), oldi = touint(oldf);
  while(oldf < f2)
  {
    __sync_bool_compare_and_swap(i, oldi, newi);
    oldf = f[0];
    oldi = touint(oldf);
  }
}

static inline void common_atomic_add(float *f, float inc)
{
  uint32_t *i = (uint32_t *)f;
  float newval;
  uint32_t newi, oldi;
  do
  {
    newval = f[0];
    oldi = touint(newval);
    newval += inc;
    newi = touint(newval);
  }
  while(!__sync_bool_compare_and_swap(i, oldi, newi));
}

static inline void common_atomic_add128(float *restrict f, const __m128 inc)
{
  unsigned __int128 *i = (unsigned __int128 *)f;
  unsigned __int128 newi, oldi;
  do
  {
    const __m128 newval = _mm_load_ps(f);
    oldi = *(unsigned __int128 *)&newval;
    _mm_add_ps(newval, inc);
    newi = *(unsigned __int128 *)&newval;
  }
  while(!__sync_bool_compare_and_swap(i, oldi, newi));
}

static inline double common_atomic_add64(double *f, double inc)
{
  uint64_t *i = (uint64_t *)f;
  double newval;
  uint64_t newi, oldi;
  do
  {
    newval = f[0];
    oldi = touint64(newval);
    newval += inc;
    newi = touint64(newval);
  }
  while(!__sync_bool_compare_and_swap(i, oldi, newi));
  return newval;
}

static inline void common_sincosf(float phi, float* sin, float* cos)
{
#ifdef __APPLE__
  *sin = sinf(phi);
  *cos = cosf(phi);
#else
  sincosf(phi, sin, cos);
#endif
}

static inline ssize_t common_readahead(int fd, off_t offset, size_t count)
{
#ifdef __APPLE__
  //??
  return 0;
#else
  return readahead(fd, offset, count);
#endif
}

static inline void common_setaffinity(int cpuid)
{
#ifdef __APPLE__
  // TODO: find out how
#else
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(cpuid, &set);
  sched_setaffinity(0, sizeof(cpu_set_t), &set);
#endif
}

#if 1 // the following block is stolen from:
/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

static inline float
common_fasterpow2 (float p)
{
  float clipp = (p < -126) ? -126.0f : p;
  union { uint32_t i; float f; } v = { (uint32_t) ( (1 << 23) * (clipp + 126.94269504f) ) };
  return v.f;
}

static inline float
common_fasterexp (float p)
{
  return common_fasterpow2 (1.442695040f * p);
}

static inline float 
common_fasterlog2 (float x)
{
  union { float f; uint32_t i; } vx = { x };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;
  return y - 126.94269504f;
}

static inline float
common_fasterlog (float x)
{
  //  return 0.69314718f * fasterlog2 (x);

  union { float f; uint32_t i; } vx = { x };
  float y = vx.i;
  y *= 8.2629582881927490e-8f;
  return y - 87.989971088f;
}
#endif
