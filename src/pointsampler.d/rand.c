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
#include "sampler.h"
#include "render.h"
#include "points.h"

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: none\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = calloc(1, sizeof(*s));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  return s;
}

float pointsampler(path_t *p, int dim)
{
  // pure random mersenne twister
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[dim];
  return points_rand(rt.points, tid);
}

void pointsampler_splat(path_t *p, mf_t value)
{
  render_splat(p, value);
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  path_init(tent, tent->index, tent->sensor.camid);
  sampler_create_path(tent);
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  path_init(tent, tent->index, tent->sensor.camid);
  path_set_pixel(tent, i, j);
  sampler_create_path(tent);
}

int pointsampler_accept(path_t *curr, path_t *tent) { return 0; }
void pointsampler_clear() {}
void pointsampler_cleanup(pointsampler_t *s)
{
  free(s->rand);
  free(s);
}
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}
void pointsampler_reset_thread(pointsampler_t *t) {}
void pointsampler_finalize(pointsampler_t *s) {}

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
