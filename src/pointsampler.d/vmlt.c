/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pointsampler.h"
#include "corona_common.h"
#include "pathspace.h"
#include "points.h"
#include "sampler.h"
#include "render.h"
#include "threads.h"
#include "pathspace/vmlt.h"

#include <stdio.h>
#include <float.h>
#include <pthread.h>
#include <stdatomic.h>

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  vmlt_t *vmlt;

  // mean image brightness estimated by all threads simultaneously
  pthread_mutex_t mutex;
  uint64_t num_samples;
  uint64_t num_zero_samples;
  double b;

  // for inverse kmlt sampling:
  fake_randoms_t *rand;
}
pointsampler_t;

#define POINTSAMPLER_INIT_SAMPLES (1ul<<20)

// makes little sense for veach mlt, we always need a single path (not a group)
void pointsampler_splat(path_t *path, float value) {}

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: metropolis light transport in path space (veach style)\n");
  fprintf(f, "           using %lu max consecutive rejects.\n", rt.pointsampler->vmlt->max_num_rejects);
  fprintf(f, "           mean image brightness %g\n", rt.pointsampler->b);
  for(int k=0;k<rt.pointsampler->vmlt->mutations;k++)
  {
    uint64_t accepted = 0, rejected = 0;
    rt.pointsampler->vmlt->mutation[k].print_info(f, rt.pointsampler->vmlt->mutation[k].data);
    for(int i=0;i<rt.num_threads;i++)
    {
      accepted += rt.pointsampler->vmlt->t[i].accepted[k];
      rejected += rt.pointsampler->vmlt->t[i].rejected[k];
    }
    fprintf(f, "           acceptance rate %2.4f%% (%lu/%lu)\n", 100.f*(float)accepted/(float)(accepted + rejected), accepted, accepted + rejected);
  }
}

void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}

void pointsampler_cleanup(pointsampler_t *s)
{
  vmlt_cleanup(s->vmlt);
  pthread_mutex_destroy(&s->mutex);
  free(s->rand);
  free(s);
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = calloc(1, sizeof(*s));
  s->vmlt = vmlt_init();
  pthread_mutex_init(&s->mutex, 0);
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  return s;
}

float pointsampler(path_t *p, const int i)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[i];
  return points_rand(rt.points, tid);
}

int pointsampler_accept(path_t *curr, path_t *tent)
{
  pointsampler_t *s = rt.pointsampler;
  vmlt_thr_t *t = s->vmlt->t + common_get_threadid();
  if((s->num_samples < POINTSAMPLER_INIT_SAMPLES) && t->large_step)
  {
    const double thr = path_throughput(tent);
    if(thr <= 0.0) // in case many zero samples come in, count them extra until num_samples with thr > 0 have been found.
      __sync_fetch_and_add(&s->num_zero_samples, 1);
    else
    {
      pthread_mutex_lock(&s->mutex);
      if(s->num_samples < POINTSAMPLER_INIT_SAMPLES)
      {
        s->num_samples++;
        if(thr == thr) // super stupid last resort NaN check
          s->b = s->b*((s->num_samples-1.0)/(double)s->num_samples) + thr/(double)s->num_samples;
        assert(s->b == s->b);
      }
      pthread_mutex_unlock(&s->mutex);
    }
  }
  return vmlt_accept(s->vmlt, curr, tent);
}

void pointsampler_finalize(pointsampler_t *s)
{
  for(int k=0;k<rt.num_threads;k++)
  {
    vmlt_thr_t *t = rt.pointsampler->vmlt->t + k;
    view_splat(t->curr_path, t->accum);
    t->accum = 0.0f;
  }
  view_set_exposure_gain(s->b*s->num_samples/(s->num_samples+s->num_zero_samples));
}

void pointsampler_clear()
{
  rt.pointsampler->b = 0;
  view_set_exposure_gain(1.0);
  rt.pointsampler->num_samples = 0;
  memset(rt.pointsampler->vmlt->t, 0, sizeof(vmlt_thr_t)*rt.num_threads);
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  pointsampler_mutate(curr, tent); // ignore pixel position
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  vmlt_mutate(rt.pointsampler->vmlt, curr, tent);

  // first time we're called, no current sample yet.
  const int tid = common_get_threadid();
  vmlt_thr_t *t = rt.pointsampler->vmlt->t + tid;
  if((t->num_mutations == 0) && (path_measurement_contribution_dwp(tent, 0, tent->length-1) > 0.)) t->acceptance = 1.0f;
  if(t->num_mutations < (1ul << 63))
    t->num_mutations++;
}

void pointsampler_enable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 1;
}

void pointsampler_disable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 0;
}

void pointsampler_set_fake_random(pointsampler_t *s, int dim, float rand)
{
  const int tid = common_get_threadid();
  s->rand[tid].rand[dim] = rand;
}

void pointsampler_prepare_frame(pointsampler_t *s) {}
void pointsampler_stop_learning(pointsampler_t *s) {}

#undef POINTSAMPLER_INIT_SAMPLES

