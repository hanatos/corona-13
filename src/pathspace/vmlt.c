#include "points.h"
#include "view.h"
#include "pathspace/vmlt.h"
#include "vmlt_registry.h"

#include <assert.h>

void vmlt_cleanup(vmlt_t *s)
{
  for(int k=0;k<s->mutations;k++)
    s->mutation[k].cleanup(s->mutation[k].data);
  free(s->t);
  free(s);
}

vmlt_t *vmlt_init()
{
  vmlt_t *s = malloc(sizeof(*s));
  memset(s, 0, sizeof(*s));
  s->t = malloc(sizeof(*s->t)*rt.num_threads);
  memset(s->t, 0, sizeof(*s->t)*rt.num_threads);
  s->mutations = 0;
  vmlt_register_all(s);
  if(s->mutations == 0)
    fprintf(stderr, "[vmlt] no mutation strategies selected! your image will be black.\n");
  // helps against stuck chains
  s->max_num_rejects = 30000;
  return s;
}

void vmlt_register(
    vmlt_t *p,
    vmlt_init_t init,
    vmlt_cleanup_t cleanup,
    vmlt_suitability_t s,
    vmlt_mutate_t m,
    vmlt_print_info_t i)
{
  vmlt_mutation_t *t = p->mutation + p->mutations;
  t->init = init;
  t->cleanup = cleanup;
  t->suitability = s;
  t->mutate = m;
  t->print_info = i;
  p->mutations++;
  t->data = t->init();
}

int vmlt_accept(vmlt_t *s, path_t *curr, path_t *tent)
{
  const int tid = common_get_threadid();
  vmlt_thr_t *t = s->t + tid;

  // compute metropolis hastings acceptance probability:
  //
  // f_tent/T(curr->tent)  /  f_curr/T(tent->curr)  =  f_tent * T(t->c) / ( f_curr * T(c->t) )
  //
  // where f is the measurement contribution and T the conditional
  // probability to sample the path.

  // reject nan acceptance (divided by zero maybe?)
  if(!(t->acceptance >= 0.0)) t->acceptance = 0.0;

  const float a = MIN(1.0f, t->acceptance);

  // not a race, 32-bit float read is atomic, the rest doesn't matter.
  const float w_tent = a;
  const float w_curr = (1.0f - a);

  if(curr->length > 0)
    t->accum += w_curr; // remember to accumulate once we jump out of it

  if((points_rand(rt.points, tid) < a || ((t->num_rejects > s->max_num_rejects) && (a > 0.0f)))
      && path_measurement_contribution_dwp(tent, 0, tent->length-1) > 0.0) // <= i hate this for performance reasons.
  { // accept
    // have to accumulate now discarded state:
    if(t->accum > 0 && curr->length > 0)
      view_splat(curr, t->accum);
    t->accum = w_tent; // sample currently under consideration
    t->accepted[t->mutation]++;

    t->num_rejects = 0;
    // make sure we don't accept a rubbish path:
    assert(tent->length > 0);
    assert(tent->lambda > 0);
    // assert(path_measurement_contribution_dwp(tent) > 0.0);
    return 1;
  }
  else
  { // reject
    // accum rejected sample, only if it's not totally bogus (avoid nan pixel coords etc):
    if(a > 0.0f) view_splat(tent, w_tent);
    t->rejected[t->mutation]++;
    t->num_rejects++;
#if 0
    // if(t->num_rejects >= s->max_num_rejects && curr->length > 3)
    if(t->num_rejects >= 2000 && !t->large_step)
    {
      fprintf(stderr, "%lu consecutive rejects a=%g (%g), mutation=%d largestep=%d pt throughputs: (%g %g)\n",
          t->num_rejects, a, t->acceptance, t->mutation, t->large_step, sampler_throughput(curr), sampler_throughput(tent));
      fprintf(stderr, "curr:\n");
      path_print(curr, stderr);
      fprintf(stderr, "tent:\n");
      path_print(tent, stderr);
    }
#endif
    return 0;
  }
}

void vmlt_mutate(vmlt_t *s, path_t *curr, path_t *tent)
{
  const int num_mutations = s->mutations;
  if(num_mutations <= 0) return;
  const int tid = common_get_threadid();
  vmlt_thr_t *t = s->t + tid;
  t->curr_path = curr;

  t->large_step = 0;

  // get suitability from all modules, sample one, call mutate on it.
  float suit[num_mutations];
  float sum = 0.0f;
  for(int k=0;k<num_mutations;k++)
  {
    sum += (suit[k] = s->mutation[k].suitability(curr, s->mutation[k].data));
    if(k) suit[k] += suit[k-1];
  }
  for(int k=0;k<num_mutations-1;k++)
    suit[k] /= sum;
  suit[num_mutations-1] = 1.001f;

  t->acceptance = 0.0f;
  const float rand = points_rand(rt.points, tid);
  for(int k=0;k<num_mutations;k++)
  {
    if(rand < suit[k])
    {
      t->acceptance = s->mutation[k].mutate(curr, tent, s->mutation[k].data);
      // mark this as a sample suitable to compute mean brightness via path_throughput()
      if(k == 0) t->large_step = 1;
      t->mutation = k; // only set after mutate() so they can access last strategy
      break;
    }
  }
}
