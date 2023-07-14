#include "corona_common.h"
#include "render.h"
#include "points.h"
#include "pointsampler.h"
#include "threads.h"
#include "pathspace.h"

typedef struct render_t
{
  atomic_int_fast32_t clear_tls; // used for tls sync
}
render_t;

typedef struct render_tls_t
{
  // storage for markov chains
  path_t path0;
  path_t path1;
  // pointers to the above, for fast swapping
  path_t *curr_path;
  path_t *tent_path;
}
render_tls_t;

static void *clear_tls(void *arg)
{
  path_init(rt_tls.render->curr_path, 0, 0);
  path_init(rt_tls.render->tent_path, 0, 0);

  // TODO: wrap this somewhat more nicely in threads.h
  // give other threads in the pool a chance to run,
  // avoid picking up another job:
  rt.render->clear_tls++;
  while(rt.render->clear_tls < rt.num_threads) sched_yield();
  return 0;
}

render_t *render_init()
{
  render_t *r = (render_t *)common_alloc(256, sizeof(render_t));
  memset(r, 0, sizeof(*r));
  return r;
}

void render_cleanup(render_t *r)
{
  free(r);
}

render_tls_t *render_tls_init()
{
  render_tls_t *r = (render_tls_t *)common_alloc(256, sizeof(render_tls_t));
  path_init(&r->path0, 0, 0);
  path_init(&r->path1, 0, 0);
  r->curr_path = &r->path0;
  r->tent_path = &r->path1;
  return r;
}

void render_tls_cleanup(render_tls_t *r)
{
  free(r);
}

void render_clear()
{
  // clear thread local storage:
  threads_t *t = rt.threads;
  rt.render->clear_tls = 0;
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, clear_tls, t);
  pthread_pool_wait(&t->pool);
}


void render_print_info(FILE *fd)
{
  fprintf(fd, "render   : global illumination\n");
}

void render_sample_path(uint64_t index)
{
  path_t *tent = rt_tls.render->tent_path;
  path_t *curr = rt_tls.render->curr_path;
  tent->index = index;
  // get reproducible results even in threaded renders (note that --batch will
  // change the results though, as anim_frame will be incremented in batches):
  points_set_state(rt.points, common_get_threadid(), index, rt.anim_frame);
  // TODO: support --tiled and call mutate_with_pixel here?
#if 0
  uint64_t w = view_width(), h = view_height();
  uint64_t frame = index / (w*h);
  uint64_t y = (index - frame * (w*h))/w;
  uint64_t x = (index - frame * (w*h) - y * w);
  pointsampler_mutate_with_pixel(curr, tent, x, y);
#else
  pointsampler_mutate(curr, tent);
#endif
  if(pointsampler_accept(curr, tent))
  {
    path_t *tmp = curr;
    rt_tls.render->curr_path = tent;
    rt_tls.render->tent_path = tmp;
  }
}

void render_splat(const path_t *p, const mf_t value)
{
  view_splat(p, value);
}
