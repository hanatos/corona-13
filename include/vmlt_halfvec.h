#pragma once

#include "points.h"
#include "pathspace/halfvec.h"

#include "shader.h"
#include "camera.h"

typedef struct vmlt_halfvec_t
{
  halfvec_stats_t *stats;
}
vmlt_halfvec_t;

void *halfvec_init()
{
  vmlt_halfvec_t *d = (vmlt_halfvec_t *)malloc(sizeof(vmlt_halfvec_t));
  memset(d, 0, sizeof(vmlt_halfvec_t));
  d->stats = (halfvec_stats_t *)malloc(sizeof(halfvec_stats_t)*rt.num_threads);
  memset(d->stats, 0, rt.num_threads*sizeof(halfvec_stats_t));
  return d;
}

static inline void halfvec_stats_accum(halfvec_stats_t *accum, const vmlt_halfvec_t *d)
{
  memset(accum, 0, sizeof(halfvec_stats_t));
  for(int k=0;k<rt.num_threads;k++)
  {
    accum->project_failed += d->stats[k].project_failed;
    accum->propagate_failed += d->stats[k].propagate_failed;
    accum->perturb_failed += d->stats[k].perturb_failed;
    accum->no_throughput += d->stats[k].no_throughput;
    accum->mutations += d->stats[k].mutations;
    accum->newton_walks += d->stats[k].newton_walks;
    accum->successful_newton_walks += d->stats[k].successful_newton_walks;
    accum->mutations_proposed += d->stats[k].mutations_proposed;
    accum->reverse_check_failed += d->stats[k].reverse_check_failed;
    accum->sum_iterations += d->stats[k].sum_iterations;
    accum->sum_successful_reduces += d->stats[k].sum_successful_reduces;
    accum->raydiff_compute_failed += d->stats[k].raydiff_compute_failed;
    accum->error_increased += d->stats[k].error_increased;
  }
}

void halfvec_cleanup(void *data)
{
  vmlt_halfvec_t *d = (vmlt_halfvec_t *)data;
  free(d->stats);
  free(d);
}

float halfvec_suitability(const path_t *p, void *data)
{
  // can't construct path from scratch, so demand initialized sample (length > 0)
  // also two vertex paths make no sense (no half vectors)
  if(p->length <= 2) return 0.0f;
  // need at least one non-specular event in between, or else no half vectors to perturb
  for(int k=1;k<p->length-1;k++) if(!(p->v[k].mode & s_specular)) return 1.0f;
  return 0.0f;
}

// float halfvec_mutate_unittest(path_t *curr, path_t *tent, void *data)
float halfvec_mutate(path_t *curr, path_t *tent, void *data)
{
  vmlt_halfvec_t *d = (vmlt_halfvec_t *)data;
  // replace that by _thread at some point
  const int tid = common_get_threadid();
  // filter stats to only view stuck chains:
  // halfvec_stats_t dummy_stats = {0};
  // halfvec_stats_t *stats = (rt.pointsampler->t[tid].num_rejects < 400) ? &dummy_stats : d->stats+tid;
  halfvec_stats_t *stats = d->stats+tid;
  stats->mutations++;

  // obtain jacobian and init half vectors on curr
  // take path curr, copy over, redistribute all vertices but the one on the eye and light
  // fix point on the sensor, light can move on a plane (we'll move it back, don't worry).

  // if curr wasn't created using a half vector mutation there's a good chance ray diffs will not be inited:
  float curr_R[9*PATHSPACE_MAX_VERTS];
  float curr_rd_u[PATHSPACE_MAX_VERTS];
  float curr_rd_v[PATHSPACE_MAX_VERTS];
  float tent_R[9*PATHSPACE_MAX_VERTS];
  float tent_rd_u[PATHSPACE_MAX_VERTS];
  float tent_rd_v[PATHSPACE_MAX_VERTS];
  // scramble volume tangent frame orientation
  curr->tangent_frame_scrambling = 0.1f + points_rand(rt.points, tid)*(0.9f-0.1f);
  for(int k=1;k<curr->length-1;k++) manifold_init(curr, k);

  double curr_dh_dx = halfvec_precache(stats, curr, 0, curr->length-1, curr_R, curr_rd_u, curr_rd_v);
  if(curr_dh_dx == 0.0) return 0.0f;

  double curr_f = halfvec_measurement(curr, 0, curr->length-1);
  if(curr_f == 0.0)
  { // rubbish path, sorry dudes.
    stats->no_throughput++;
    return 0.0f;
  }

  // create a new path by perturbing half vectors:
  const float vol_pdf = halfvec_perturb(stats, curr, tent, 0, curr->length-1, curr_R, curr_rd_u, curr_rd_v);
  if(vol_pdf <= 0.0f)
  {
    stats->perturb_failed++;
    return 0.0f;
  }

  // check if we can actually walk back, to satisfy detailed balance
  const float vol_rpdf = halfvec_reverse_check(stats, curr, tent, 0, curr->length-1);
  if(vol_rpdf <= 0.0f) return 0.0f;

  // compute transition probabilities:
  const double p_tent = halfvec_pdf_perturb(stats, curr, tent, 0, curr->length-1, curr_R, curr_rd_u, curr_rd_v);
  if(!(p_tent > 0.0f))
  {
    stats->no_throughput++;
    return 0.0f;
  }

  // computing the reversre pdf requires ray differentials to be initialized on tent
  double tent_dh_dx = halfvec_compute_raydifferentials(stats, tent, 0, curr->length-1, tent_R, tent_rd_u, tent_rd_v);
  if(tent_dh_dx == 0.0) return 0.0f;
  const double p_curr = halfvec_pdf_perturb(stats, tent, curr, 0, curr->length-1, tent_R, tent_rd_u, tent_rd_v);
  if(!(p_curr > 0.0f))
  {
    stats->reverse_check_failed++;
    return 0.0f;
  }


  // compute acceptance as measurement contribution times transform from halfvector space
  // to vertex area measure (multiply the jacobian).
  double tent_f = halfvec_measurement(tent, 0, tent->length-1);
  if(tent_f <= 0.0f)
    stats->no_throughput++;
  stats->mutations_proposed++;

  // compute acceptance as f->/T->  /  f<-/T<-
  return tent_f/p_tent / (curr_f / p_curr);
}

void halfvec_print_info(FILE *f, void *data)
{
  vmlt_halfvec_t *d = (vmlt_halfvec_t *)data;
  halfvec_stats_t stats;
  halfvec_stats_accum(&stats, d);
  fprintf(f, "         : halfvector space mutation\n");
  fprintf(f, "           failed project         %.02f%%\n", 100.0f*stats.project_failed/(float)stats.mutations);
  fprintf(f, "           failed propagate       %.02f%%\n", 100.0f*stats.propagate_failed/(float)stats.mutations);
  fprintf(f, "           no throughput          %.02f%%\n", 100.0f*stats.no_throughput/(float)stats.mutations);
  fprintf(f, "           failed perturbation    %.02f%%\n", 100.0f*stats.perturb_failed/(float)stats.mutations);
  fprintf(f, "             did not converge     %.02f%%\n", 100.0f*(stats.newton_walks-stats.successful_newton_walks)/(float)stats.newton_walks);
  fprintf(f, "             average iterations   %.02f\n", stats.sum_iterations/(float)stats.newton_walks);
  fprintf(f, "             successful step red. %.02f%%\n", 100.0f*stats.sum_successful_reduces/(float)stats.sum_iterations);
  fprintf(f, "             error increased      %.02f%%\n", 100.0f*stats.error_increased/(float)stats.sum_iterations);
  fprintf(f, "           non-reversible         %.02f%%\n", 100.0f*stats.reverse_check_failed/(float)stats.mutations);
  fprintf(f, "           raydiff compute failed %.02f%%\n", 100.0f*stats.raydiff_compute_failed/(float)stats.mutations);
  fprintf(f, "           successful mutations   %.02f%%\n", 100.0f*stats.mutations_proposed/(float)stats.mutations);
}
