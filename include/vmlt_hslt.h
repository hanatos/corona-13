#pragma once

#include "points.h"
#include "pathspace/hslt.h"

#include "shader.h"
#include "camera.h"

typedef struct vmlt_hslt_t
{
  halfvec_stats_t *stats;
}
vmlt_hslt_t;

void *hslt_init()
{
  vmlt_hslt_t *d = (vmlt_hslt_t *)malloc(sizeof(vmlt_hslt_t));
  memset(d, 0, sizeof(vmlt_hslt_t));
  d->stats = (halfvec_stats_t *)malloc(sizeof(halfvec_stats_t)*rt.num_threads);
  memset(d->stats, 0, rt.num_threads*sizeof(halfvec_stats_t));
  return d;
}

static inline void hslt_stats_accum(halfvec_stats_t *accum, const vmlt_hslt_t *d)
{
  memset(accum, 0, sizeof(halfvec_stats_t));
  for(int k=0;k<rt.num_threads;k++)
  {
    accum->newton_walks += d->stats[k].newton_walks;
    accum->successful_newton_walks += d->stats[k].successful_newton_walks;
    accum->project_failed += d->stats[k].project_failed;
    accum->propagate_failed += d->stats[k].propagate_failed;
    accum->perturb_failed += d->stats[k].perturb_failed;
    accum->no_throughput += d->stats[k].no_throughput;
    accum->mutations += d->stats[k].mutations;
    accum->mutations_proposed += d->stats[k].mutations_proposed;
    accum->reverse_check_failed += d->stats[k].reverse_check_failed;
    accum->sum_iterations += d->stats[k].sum_iterations;
    accum->sum_successful_reduces += d->stats[k].sum_successful_reduces;
    accum->raydiff_compute_failed += d->stats[k].raydiff_compute_failed;
    accum->error_increased += d->stats[k].error_increased;
  }
}

void hslt_cleanup(void *data)
{
  vmlt_hslt_t *d = (vmlt_hslt_t *)data;
  free(d->stats);
  free(d);
}

float hslt_suitability(const path_t *p, void *data)
{
  // can't construct path from scratch, so demand initialized sample (length > 0)
  if(p->length < 2) return 0.0f;
  return 1.0f;
}

float hslt_mutate(
    path_t *curr,       // constant current sample
    path_t *tent,       // tentative sample
    void *data)
{
  vmlt_hslt_t *d = (vmlt_hslt_t *)data;
  // replace that by _thread at some point
  const int tid = common_get_threadid();
  d->stats[tid].mutations++;

  // scramble volume tangent frame orientation
  curr->tangent_frame_scrambling = 0.1f + points_rand(rt.points, tid)*(0.9f-0.1f);
  for(int k=1;k<curr->length-1;k++) manifold_init(curr, k);

  // we need to init half vectors to be able to perturb them!
  // init current path's half vectors, and constraint derivatives.
  curr->v[0].diffgeo.type = s_pinned_position;
  for(int i=1;i<curr->length;i++)
    curr->v[i].diffgeo.type = s_free;
  if(manifold_compute_tangents(curr, 0, curr->length-1) == 0.0)
    return 0.0f;

  // first very wastefully copy the whole path
  // TODO: really not necessary, could only do vertex by vertex (multichain inits them anyhow)
  // TODO: maybe do this, but then construct multichain on different, empty path
  // and use this to converge the halfvector space part in the same index range?
  *tent = *curr;

  // determine vertex indices a=0, b=?, c=?:

  int b = tent->length - 1; // default to all lens perturbation
  float pdf_b_tent = 1.0f;
  float b_cdf[PATHSPACE_MAX_VERTS];
  if(tent->length > 2)
  {
    hslt_get_b_cdf(tent, 0, b_cdf);
    b = sample_cdf(b_cdf, tent->length, points_rand(rt.points, tid));
    pdf_b_tent = b ? b_cdf[b] - b_cdf[b-1] : b_cdf[0];
  }

  float c_cdf[PATHSPACE_MAX_VERTS];
  float pdf_c_tent = 1.0f;
  int c = tent->length - 1;
  if(b < tent->length - 1)
  {
    hslt_get_c_cdf(tent, b, c_cdf);
    c = b + 1 + sample_cdf(c_cdf+b+1, tent->length-b-1, points_rand(rt.points, tid));
    pdf_c_tent = c ? c_cdf[c] - c_cdf[c-1] : c_cdf[0];
  }

  assert(!(tent->v[b].mode & s_specular));
  assert(!(tent->v[c].mode & s_specular));


  // track pdf ratio T(t->c)/T(c->t) in vertex area measure (used for camera vertex)
  if(b > 0)
  {
    float g1, g2;
    const float r1 = points_rand(rt.points, tid);
    const float r2 = points_rand(rt.points, tid);
    sample_gaussian(r1, r2, &g1, &g2);
    // step how many pixels?
    // sub-manifold boundaries here:
    const float stepsize = HALFVEC_MUTATION_STEP;
    tent->sensor.pixel_i += stepsize * g1;
    tent->sensor.pixel_j += stepsize * g2;
    // mutate point on aperture, after mutating outgoing direction (pixel)
    tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.01f);
    tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.01f);
    tent->sensor.pixel_set = 1;
    tent->sensor.aperture_set = 1;
  }

  double T = hslt_perturb(curr, tent, b, c, 0, d->stats+tid);
  if(T <= 0.0f)
  {
    d->stats[tid].perturb_failed++;
    return 0.0f;
  }

  const int to = tent->length - curr->length;
  if(tent->length > 2)
  {
    hslt_get_b_cdf(tent, 0, b_cdf);
    const float pdf_b_curr = (b+to) ? b_cdf[b+to] - b_cdf[b+to-1] : b_cdf[0];
    hslt_get_c_cdf(tent, b+to, c_cdf);
    const float pdf_c_curr = c+to ? c_cdf[c+to] - c_cdf[c+to-1] : c_cdf[0];

    T *= pdf_b_curr * pdf_c_curr / (pdf_b_tent * pdf_c_tent);
  }

  // return acceptance as f->/T->  /  f<-/T<-
  return T;
}

void hslt_print_info(FILE *f, void *data)
{
  vmlt_hslt_t *d = (vmlt_hslt_t *)data;
  halfvec_stats_t stats;
  hslt_stats_accum(&stats, d);
  fprintf(f, "         : extended half vector space mutation\n");
  fprintf(f, "           newtonian walks        %.2f%%\n", 100.0f*stats.newton_walks/(float)stats.mutations);
  fprintf(f, "           failed project         %.2f%%\n", 100.0f*stats.project_failed/(float)stats.mutations);
  // fprintf(f, "           failed propagate       %.02f%%\n", 100.0f*stats.propagate_failed/(float)stats.mutations);
  fprintf(f, "           no throughput          %.02f%%\n", 100.0f*stats.no_throughput/(float)stats.mutations);
  fprintf(f, "           failed perturbation    %.02f%%\n", 100.0f*stats.perturb_failed/(float)stats.mutations);
  fprintf(f, "             did not converge     %.02f%%\n", 100.0f*(stats.newton_walks-stats.successful_newton_walks)/(float)stats.newton_walks);
  fprintf(f, "             average iterations   %.02f\n", stats.sum_iterations/(float)stats.successful_newton_walks);
  fprintf(f, "             successful step red. %.02f%%\n", 100.0f*stats.sum_successful_reduces/(float)stats.sum_iterations);
  fprintf(f, "             error increased      %.02f%%\n", 100.0f*stats.error_increased/(float)stats.sum_iterations);
  fprintf(f, "           non-reversible         %.02f%%\n", 100.0f*stats.reverse_check_failed/(float)stats.newton_walks);
  fprintf(f, "           raydiff compute failed %.02f%%\n", 100.0f*stats.raydiff_compute_failed/(float)stats.newton_walks);
  fprintf(f, "           successful mutations   %.02f%%\n", 100.0f*stats.mutations_proposed/(float)stats.mutations);
}
