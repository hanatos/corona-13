#include "corona_common.h"
#include "render.h"
#include "points.h"
#include "pointsampler.h"
#include "threads.h"
#include "pathspace.h"
#include "pathspace/vmlt.h"
#include "pathspace/tech.h"
#include "sampler.h"
#include "filter.h"

typedef struct render_t
{
  atomic_int_fast32_t clear_tls; // used for tls sync

  vmlt_t *vmlt;

  // cascades of throughput buffers:
  int stops;           // how many stops difference between buffers
  int num_buffers;     // how many buffers
  uint64_t buf_width;  // dimensions of buffer
  uint64_t buf_height;
  float **throughputs; // the buffers themselves, 1 channel per pixel, half res
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
  render_t *r = malloc(sizeof(*r));
  memset(r, 0, sizeof(*r));
  r->stops = 1;
  r->num_buffers = 8;
  r->buf_width  = view_width()/2;
  r->buf_height = view_height()/2;
  r->throughputs = malloc(sizeof(float*)*r->num_buffers);
  for(int k=0;k<r->num_buffers;k++)
    r->throughputs[k] = calloc(r->buf_width*r->buf_height, sizeof(float));
  r->vmlt = vmlt_init();
  return r;
}

void render_cleanup(render_t *r)
{
#if 1 // XXX DEBUG
  for(int k=0;k<r->num_buffers;k++)
  {
    char filename[1024];
    snprintf(filename, sizeof(filename), "debug_%02d.pfm", k);
    FILE *f = fopen(filename, "wb");
    fprintf(f, "PF\n%lu %lu\n-1.0\n", r->buf_width, r->buf_height);
    for(uint64_t i=0;i<r->buf_width*r->buf_height;i++)
      for(int j=0;j<3;j++)
        fwrite(r->throughputs[k]+i, 1, sizeof(float), f);
    fclose(f);
  }
#endif 
  vmlt_cleanup(r->vmlt);
  free(r);
}

render_tls_t *render_tls_init()
{
  render_tls_t *r = (render_tls_t *)malloc(sizeof(render_tls_t));
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
  fprintf(fd, "render   : erpt\n");
  // TODO: dump stats about vmlt
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

static inline int _get_mutations(const float v, int *buf)
{
  const float logval = logf(v+1.0f)/logf(2.0f)/2.0f;
  const int b = CLAMP((int)logval, 0, rt.render->num_buffers-1);
  if(buf) *buf = b;
  // avoid super long run times for first frames when throughput
  // buffers are still empty:
  if(rt.frames < 2) return 1<<b;
  // somehow this also needs to compensate for autocorrelation of
  // the markov chain:
  return 8<<b;
}

void render_splat(const path_t *p, const float value)
{
  if(!(value > 0)) return;
  if(p->length < 3)
  {
    view_splat(p, value);
    return;
  }
  // splat 1 into throughput buffer cascade selected by value
  // if readout of buffer < trust value:
  //  run erpt chain to bring down brightness into, uhm.. first cascade?
  int buf = 0;
  int mutations = _get_mutations(value, &buf);
  // TODO: splat fractional parts in both adjacent buffers?
  float col = 1.0f;
  const float x = p->sensor.pixel_i/2.0f, y = p->sensor.pixel_j/2.0f;
  // TODO: really need a pixel filter?
  filter_splat(x, y, &col,
      rt.render->throughputs[buf], 1, 0, 1,
      rt.render->buf_width, rt.render->buf_height);

  // TODO: more clever readout? bilinear at least? look at bilateral grid i guess
  int i = CLAMP((int)x, 0, rt.render->buf_width-1);
  int j = CLAMP((int)y, 0, rt.render->buf_height-1);
  const float trust = rt.render->throughputs[buf][i + j*rt.render->buf_width];

  if(trust < 2.0f)
  {
    path_t path1, path2;
    path_t *curr = &path1, *tent = &path2;

    // TODO: if from the light, reverse path
    // TODO: have path_copy() that only copies the active vertices/edges:
    *curr = *p;

    // fprintf(stderr, "mutations %d trust %g\n", mutations, trust);
    float accum = 0.0f;
    // mutate
    for(int m=0;m<mutations;m++)
    { // stolen metropolis acceptance logic from pointsampler_vmlt.c:
      vmlt_mutate(rt.render->vmlt, curr, tent);
      const int tid = common_get_threadid();
      vmlt_thr_t *t = rt.render->vmlt->t + tid;
      if(!(t->acceptance >= 0.0)) t->acceptance = 0.0;

      float a = MIN(1.0f, t->acceptance);
      if(a > 0.0f)
      {
        // weight: their buffer #mutations / mutations
        const float tent_throughput = sampler_mis_weight(tent) *
            path_measurement_contribution_dx(tent, 0, tent->length-1)/
            path_tech_pdf_as_sampled(tent);
        const int mutations_rev = _get_mutations(tent_throughput, 0);
        // fprintf(stderr, "mut %d %d\n", mutations, mutations_rev);
        // reject if this falls out of our class:
        if(mutations != mutations_rev) a = 0.0f; // TODO: weighted/probabilistic?
        const float w = mutations_rev / (float)(mutations*mutations);
        const float w_tent = w * value * a;
        const float w_curr = w * value * (1.0f - a);
        if(curr->length > 0) accum += w_curr;
        if(points_rand(rt.points, common_get_threadid()) < a)
        { // accept
          if(accum > 0.0 && curr->length > 0)
            view_splat(curr, accum);
          accum = w_tent;
          // swap paths:
          path_t *tmp = curr;
          curr = tent;
          tent = tmp;
        }
        else
        { // reject
          view_splat(tent, w_tent);
        }
      }
    }
    // accumulate the rest:
    if(accum > 0.0 && curr->length > 0)
      view_splat(curr, accum);
  }
  else // got good trust in the sample:
  {
    // go ahead and splat as usual
    view_splat(p, value);
  }
}

