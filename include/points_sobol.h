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
#ifndef POINTS_H
#define POINTS_H


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#ifdef __SSE2__
  #include <emmintrin.h>
#endif

extern double drand48(void);

typedef struct sobol_t
{
  unsigned int *matrices;
  unsigned int nDimensions;
}
sobol_t;

typedef struct sobol_state_t
{
  unsigned int *state;
  unsigned int index;
  unsigned int c;
}
sobol_state_t;

static inline void sobol_init(sobol_t *s, const unsigned int nDimensions, const char * const matrixDataFilename)
{
  s->matrices = NULL;
  s->nDimensions = nDimensions;
  if (nDimensions & 3)
  {
    fprintf(stderr, "nDimensions must be a multiple of four!\n");
    return;
  }

  gzFile f;
  if(matrixDataFilename == NULL) f = gzopen("data/points_sobol.gz", "rb");
  else f = gzopen(matrixDataFilename, "rb");
  if (!f)
  {
    fprintf(stderr, "Error opening %s !\n", matrixDataFilename);
    return;
  }

  unsigned int maxDim;
  gzread(f, &maxDim, sizeof(unsigned int));

  if (nDimensions > maxDim)
  {
    fprintf(stderr, "%s does not supply %d dimensions!\n", matrixDataFilename, nDimensions);
    gzclose(f);
    return;
  }

  s->matrices = (unsigned int *)common_alloc(64, sizeof(unsigned int) * s->nDimensions * 32);

  for (unsigned int k = 0; k < 32; k++)
  {
    gzseek(f, (1 + k * maxDim) * sizeof(unsigned int), SEEK_SET);
    gzread(f, s->matrices + k * s->nDimensions, sizeof(unsigned int) * s->nDimensions);
  }

  gzclose(f);
}

static inline void sobol_set_state(const sobol_t * const m, sobol_state_t * s, unsigned int pointIndex)
{
  s->index = pointIndex;
  s->c = 0;
  pointIndex ^= pointIndex >> 1; // Gray code representation

  memset(s->state, 0, sizeof(unsigned int) * m->nDimensions);

  const unsigned int *matrix = m->matrices;
  for (; pointIndex; pointIndex >>= 1, matrix += m->nDimensions)
    if (pointIndex & 1)
      for (unsigned int d = 0; d < m->nDimensions; d++)
        s->state[d] ^= matrix[d];
}
    
static inline void sobol_init_state(sobol_t *m, sobol_state_t *s, const unsigned int index)
{
  s->state = NULL;
  s->index = 0;
  s->c = 0;
  s->state = (unsigned int *)common_alloc(64, sizeof(unsigned int) * (m->nDimensions));
  sobol_set_state(m, s, index);
}

static inline void sobol_cleanup_state(sobol_state_t *s)
{
  if(s->state) free(s->state);
  s->state = 0;
}

static inline void sobol_cleanup(sobol_t *s)
{
  if(s->matrices) free(s->matrices);
  s->matrices = 0;
}

static inline void sobol_next(const sobol_t * const m, sobol_state_t *s)
{
  // next point
  if(!++s->index)
  {
    sobol_set_state(m, s, 0);
    return;
  }

  // find the bit that has changed (Gray code)
#ifdef __ICC
  const int pos = _bit_scan_forward(s->index);
#else
#if GCC_VERSION >= 3004
  const int pos = __builtin_ctz(s->index);
#else
  int pos = 0;
  for (int mask = 1; !(s->index & mask); ++pos, mask <<= 1);
#endif
#endif

  // update the state
#ifdef __SSE2__
  __m128i *m2 = (__m128i *) (m->matrices + pos * m->nDimensions);
  __m128i *s2 = (__m128i *) (s->state);
  const unsigned int maxD = m->nDimensions >> 2;
  for (unsigned int d = 0; d < maxD; d++)
    s2[d] = _mm_xor_si128(s2[d], m2[d]);
#else
  unsigned int *t = m->matrices + pos * m->nDimensions;
  for (unsigned int d = 0; d < m->nDimensions; d++)
    s->state[d] ^= t[d];
#endif
}

static inline double sobol_component(const sobol_state_t * const s, const unsigned int c)
{
  const unsigned long long dreggn = 1ULL << 32;
  const double scale = 1. / dreggn;
  return s->state[c] * scale;
}

typedef struct points_t
{
  sobol_t m;
  sobol_state_t **s;
  unsigned int num;
}
points_t;

static inline void points_init(points_t *p, const unsigned int num_threads)
{
  const int dim = 20; // max number of dimensions before going random
  p->s = (sobol_state_t **)malloc(sizeof(sobol_state_t *)*num_threads);
  p->s[0] = (sobol_state_t *)malloc(sizeof(sobol_state_t)*num_threads);
  p->num = num_threads;
  sobol_init(&(p->m), dim, NULL);
  sobol_init_state(&(p->m), p->s[0], 0);
  for(unsigned int i=1;i<num_threads;i++)
  {
    p->s[i] = p->s[i-1] + 1;
    sobol_init_state(&(p->m), p->s[i], 0);
  }
}

static inline void points_cleanup(points_t *p)
{
  sobol_cleanup(&(p->m));
  for(unsigned int i=0;i<p->num;i++)
  {
    sobol_cleanup_state(p->s[i]);
  }
  free(p->s[0]);
  free(p->s);
}

static inline void points_set_state(points_t *p, const unsigned int thread_num, const unsigned int i)
{
  sobol_set_state(&(p->m), p->s[thread_num], i);
}

static inline void points_next(points_t *p, const unsigned int thread_num)
{
  sobol_next(&(p->m), p->s[thread_num]);
  p->s[thread_num]->c = 0;
}

static inline double points_rand(points_t *p, const unsigned int thread_num)
{
  if(++ p->s[thread_num]->c >= p->m.nDimensions) return drand48();
  return sobol_component(p->s[thread_num], p->s[thread_num]->c);
}

static inline double points_get(points_t *p, const unsigned int thread_num, const unsigned int c)
{
  //assert(c < p->m.nDimensions);
  if(c >= p->m.nDimensions) return drand48();
  return sobol_component(p->s[thread_num], c);
}

static inline void points_print_info(FILE *fd)
{
  fprintf(fd, "points   : sobol sequence\n");
}

#endif
