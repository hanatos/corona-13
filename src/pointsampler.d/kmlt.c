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
#include "sampler_common.h"
#include "view.h"
#include "corona_common.h"
#include "pathspace.h"
#include "points.h"
#include <stdio.h>
#include <float.h>

// some parameters:
// samples to estimate b
#define POINTSAMPLER_INIT_SAMPLES (1ul<<20)
// probability for a large step mutation
#define POINTSAMPLER_P_LARGE_STEP 0.2f
#define POINTSAMPLER_MAX_RANDS_PER_VERT 10
// leave space for bdpt (two paths)
#define POINTSAMPLER_NUM_DIMS 2*(PATHSPACE_MAX_VERTS*POINTSAMPLER_MAX_RANDS_PER_VERT)

typedef struct pointsampler_contribution_t
{
  float i,j;
  int path_length;
  int camid;
  mf_t lambda;
  mf_t throughput, throughput_mapped;
}
pointsampler_contribution_t;

typedef struct pointsampler_thr_t
{
  int mutated_samples;
  float temperature, mean_temperature;
  mf_t curr_throughput, tent_throughput, accum_throughput;
  mf_t curr_throughput_mapped, tent_throughput_mapped;
  int64_t num_samples, num_rejects, num_accepts, num_mutations;
  double b, b_mapped;
  int curr_num_contribs, tent_num_contribs;
  // two buffers and bdpt potentially creates #V^2 contributions.
  struct pointsampler_contribution_t contrib_buf[2*PATHSPACE_MAX_VERTS*PATHSPACE_MAX_VERTS];
  struct pointsampler_contribution_t *curr_contrib, *tent_contrib;
  int large_step;
  float p_large_step;
  float *rand_buf;     // buffers to be swapped.
  float *tent_rand;    // tentative state
  float *curr_rand;    // current (old) sample
}
pointsampler_thr_t;

typedef struct pointsampler_t
{
  pointsampler_thr_t *t;
}
pointsampler_t;


void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: metropolis as in kelemen and szirmay-kalos\n");
  fprintf(f, "           large step probability %f\n", POINTSAMPLER_P_LARGE_STEP);
  double mean_temp = 0.0;
  double acceptance = 0.0;
  double b = 0.0;
  for(int k=0;k<rt.num_threads;k++)
  {
    mean_temp += rt.pointsampler->t[k].mean_temperature / (double)rt.pointsampler->t[k].num_mutations;
    acceptance += rt.pointsampler->t[k].num_accepts / (double)rt.pointsampler->t[k].num_mutations;
    b += rt.pointsampler->t[k].b / (double)rt.num_threads;
  }
  mean_temp /= rt.num_threads;
  acceptance /= rt.num_threads;
  // fprintf(f, "           mean temperature %f\n", mean_temp);
  fprintf(f, "           mean acceptance %.02f%%\n", 100.0*acceptance);
  fprintf(f, "           mean image brightness %g\n", b);
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->t = (pointsampler_thr_t *)malloc(sizeof(pointsampler_thr_t)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
  {
    s->t[k].rand_buf = (float *)malloc(2*sizeof(float)*POINTSAMPLER_NUM_DIMS);
    memset(s->t[k].rand_buf, 0, 2*sizeof(float)*POINTSAMPLER_NUM_DIMS);
    s->t[k].tent_rand = s->t[k].rand_buf;
    s->t[k].curr_rand = s->t[k].rand_buf + POINTSAMPLER_NUM_DIMS;
    s->t[k].curr_rand[s_dim_lambda] = -0.666f;
    s->t[k].tent_throughput = mf_set1(0.0f);
    s->t[k].accum_throughput = mf_set1(0.0f);
    s->t[k].tent_throughput_mapped = mf_set1(0.0f);
    s->t[k].temperature = 0.0f;
    s->t[k].mean_temperature = 0.0f;
    s->t[k].curr_throughput = mf_set1(FLT_MIN);
    s->t[k].curr_throughput_mapped = mf_set1(FLT_MIN);
    s->t[k].curr_contrib = s->t[k].contrib_buf;
    s->t[k].tent_contrib = s->t[k].contrib_buf + 2*PATHSPACE_MAX_VERTS;
    s->t[k].curr_num_contribs = 0;
    s->t[k].tent_num_contribs = 0;
    s->t[k].large_step = 0;
    s->t[k].num_samples = 0;
    s->t[k].b = mf_set1(0.0f);
    s->t[k].b_mapped = mf_set1(0.0f);
    s->t[k].num_rejects = 0;
    s->t[k].p_large_step = POINTSAMPLER_P_LARGE_STEP;
  }
  return s;
}

void pointsampler_set_large_step(pointsampler_t *t, float p_large_step)
{
  const int tid = common_get_threadid();
  t->t[tid].p_large_step = p_large_step;
}

void pointsampler_finalize(pointsampler_t *s)
{
  // last sample:
  path_t path;
  path_init(&path, 0, 0);
  float b = 0.0f;
  for(int i=0;i<rt.num_threads;i++)
  {
    pointsampler_thr_t *t = s->t + i;
    if(t->accum_throughput > 0.0f) for(int k=0;k<t->curr_num_contribs;k++)
    {
      path.sensor.pixel_i = t->curr_contrib[k].i;
      path.sensor.pixel_j = t->curr_contrib[k].j;
      path.lambda = t->curr_contrib[k].lambda;
      path.length = t->curr_contrib[k].path_length;
      path.sensor.camid = t->curr_contrib[k].camid;
      view_splat(&path, t->accum_throughput * t->curr_contrib[k].throughput/t->curr_throughput);
    }
    t->accum_throughput = 0.0;
    b += t->b / rt.num_threads;
  }
  view_set_exposure_gain(b);
}

void pointsampler_clear()
{
  for(int k=0;k<rt.num_threads;k++)
  {
    rt.pointsampler->t[k].num_samples = 0;
    rt.pointsampler->t[k].num_accepts = 0;
    rt.pointsampler->t[k].num_mutations = 0;
    rt.pointsampler->t[k].temperature = 0.f;
    rt.pointsampler->t[k].mean_temperature = 0.f;
    rt.pointsampler->t[k].b = 0.0f;
    rt.pointsampler->t[k].b_mapped = 0.0f;
    rt.pointsampler->t[k].tent_throughput = 0.0f;
    rt.pointsampler->t[k].accum_throughput = 0.0f;
    rt.pointsampler->t[k].curr_throughput = FLT_MIN;
    rt.pointsampler->t[k].tent_throughput_mapped = 0.0f;
    rt.pointsampler->t[k].curr_throughput_mapped = FLT_MIN;
    rt.pointsampler->t[k].num_rejects = 0;
  }
  view_set_exposure_gain(1.0);
}

void pointsampler_cleanup(pointsampler_t *s)
{
  for(int k=0;k<rt.num_threads;k++) free(s->t[k].rand_buf);
  free(s->t);
  free(s);
}

float pointsampler(path_t *p, const int i)
{
  int v = p->length;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = rt.pointsampler->t + tid;
  const int cnt = POINTSAMPLER_MAX_RANDS_PER_VERT;
  const int end = p->v[v].rand_beg;
  if(end + i >= t->mutated_samples)
  {
    const int old = t->mutated_samples;
    t->mutated_samples = end + 3*cnt;
    if(t->mutated_samples > PATHSPACE_MAX_VERTS*cnt) t->mutated_samples = PATHSPACE_MAX_VERTS*cnt;
    for(int k=old;k<t->mutated_samples;k++)
      t->tent_rand[k] = points_rand(rt.points, tid);
  }
  assert(end + i < POINTSAMPLER_NUM_DIMS);
  return t->tent_rand[end + i];
}

static inline float tonemap(float x)
{
#if 0
  // log space, expect mean image brightness somewhere around [0,1] with is taken into account
  if(x <= 0.0) return 0.0f;
  // XXX const float tm = common_fasterlog2(1.0f + rt.cam->iso/100.0f * x);
  const float tm = common_fasterlog2(1.0f + x);
  assert(tm >= 0.0);
  return tm;
#endif

#if 0
  // return pdf of imaginary light tracing path (intersecting aperture by chance)
  double pdf = 1.0;
  if(path->v[0].mode & s_sensor)
    for(int k=0;k<path->length;k++)
      pdf *= path_pdf_extend_adjoint(path, k);
  else if(path->v[0].mode & s_emit)
    for(int k=0;k<path->length;k++)
      pdf *= path_pdf_extend(path, k);
  return x * pdf/path_pdf(path);
#endif

  // off:
  return x;
}

void pointsampler_splat(path_t *path, float value)
{ // wrapper for view_splat, postponing tentative sample.
  if(!(value > 0.0 && value < FLT_MAX)) return;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = rt.pointsampler->t + tid;
  t->tent_throughput += value;
  t->tent_throughput_mapped += tonemap(value);
  const int n = t->tent_num_contribs;
  t->tent_contrib[n].i = path->sensor.pixel_i;
  t->tent_contrib[n].j = path->sensor.pixel_j;
  t->tent_contrib[n].lambda = path->lambda;
  t->tent_contrib[n].path_length = path->length;
  t->tent_contrib[n].camid = path->sensor.camid;
  t->tent_contrib[n].throughput = value;
  t->tent_contrib[n].throughput_mapped = tonemap(value);
  t->tent_num_contribs++;
}

int pointsampler_accept(path_t *curr, path_t *tent)
{
  pointsampler_t *s = rt.pointsampler;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = s->t + tid;
  t->num_rejects++;
  t->num_mutations++;

  path_t path;
  path_init(&path, 0, 0);
  const float a = fminf(1.0f, t->tent_throughput_mapped/t->curr_throughput_mapped);
  if((t->num_samples < POINTSAMPLER_INIT_SAMPLES) && t->large_step)
  {
    t->num_samples++;
    t->b = t->b*((t->num_samples-1.0)/(double)t->num_samples) + t->tent_throughput/(double)t->num_samples;
    t->b_mapped = t->b_mapped*((t->num_samples-1.0)/(double)t->num_samples) + t->tent_throughput_mapped/(double)t->num_samples;
  }
  // divide out the wrong pdf (mapped throughput / mapped mean image brightness) and multiply with the correct but sub-optimal one
  const float w_tent = a         ;// * t->tent_throughput * t->b_mapped / (t->tent_throughput_mapped * t->b);
  const float w_curr = (1.0f - a);// * t->curr_throughput * t->b_mapped / (t->curr_throughput_mapped * t->b);

  t->accum_throughput += w_curr; // remember to accumulate once we jump out of it

  if((points_rand(rt.points, tid) < a) || (t->num_rejects > 40000 && a > 0.0))
  { // accept
    t->temperature = MAX(0, t->temperature-1);
    // have to accumulate now discarded state:
    if(t->curr_throughput > FLT_MIN)
    {
      for(int k=0;k<t->curr_num_contribs;k++)
      {
        path.sensor.pixel_i = t->curr_contrib[k].i;
        path.sensor.pixel_j = t->curr_contrib[k].j;
        path.lambda = t->curr_contrib[k].lambda;
        path.length = t->curr_contrib[k].path_length;
        path.sensor.camid = t->curr_contrib[k].camid;
        view_splat(&path, t->accum_throughput * t->curr_contrib[k].throughput/t->curr_throughput);
      }
    }
    t->accum_throughput = w_tent; // one sample currently under consideration

    float *tmp = t->tent_rand;
    t->tent_rand = t->curr_rand;
    t->curr_rand = tmp;
    pointsampler_contribution_t *tmpc = t->tent_contrib;
    t->tent_contrib = t->curr_contrib;
    t->curr_contrib = tmpc;
    t->curr_throughput = t->tent_throughput;
    t->curr_throughput_mapped = t->tent_throughput_mapped;
    t->curr_num_contribs = t->tent_num_contribs;
    t->num_rejects = 0;
    t->num_accepts++;
    return 1;
  }
  else
  { // reject
    if(t->num_rejects == 10) t->temperature = 100;
    // accum rejected sample
    for(int k=0;k<t->tent_num_contribs;k++)
    {
      path.sensor.pixel_i = t->tent_contrib[k].i;
      path.sensor.pixel_j = t->tent_contrib[k].j;
      path.lambda = t->tent_contrib[k].lambda;
      view_splat(&path, w_tent * t->tent_contrib[k].throughput/t->tent_throughput);
    }
    return 0;
  }
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  pointsampler_mutate_with_pixel(curr, tent, -1, -1);
}

void pointsampler_mutate_with_pixel(path_t *curr_path, path_t *tent_path, float x, float y)
{
  pointsampler_t *s = rt.pointsampler;
  const int tid = common_get_threadid();
  pointsampler_thr_t *t = s->t + tid;
  t->large_step = (t->curr_throughput <= FLT_MIN) || ((points_rand(rt.points, tid) < t->p_large_step) ? 1 : 0);

  if(t->large_step)
  {
    t->mutated_samples = MIN(50, POINTSAMPLER_MAX_RANDS_PER_VERT*PATHSPACE_MAX_VERTS);
    for(int k=0;k<t->mutated_samples;k++)
      t->tent_rand[k] = points_rand(rt.points, tid);
  }
  else
  {
#if 0
    for(int k=0;k<curr_path->length;k++)
    {
      // fprintf(stderr, "v[%d] beg %d cnt %d\n", k, curr_path->v[k].rand_beg, curr_path->v[k].rand_cnt);
      for(int i=k*10;i<10*(k+1);i++)
      // for(int i=curr_path->v[k].rand_beg;i<curr_path->v[k].rand_beg + curr_path->v[k].rand_cnt;i++)
    // for(int i=0;i<POINTSAMPLER_MAX_RANDS_PER_VERT*PATHSPACE_MAX_VERTS;i++)
        t->rand[i] = sample_mutate_rand(t->curr_rand[i], points_rand(rt.points, tid), 0.01f);
    }
#else
    assert(t->num_samples > 0);
    assert(t->curr_throughput > FLT_MIN);
    assert(t->curr_rand[s_dim_lambda] > -0.666f);
    t->tent_rand[s_dim_time]   = sample_mutate_rand(t->curr_rand[s_dim_time],    points_rand(rt.points, tid), 0.1f);
    t->tent_rand[s_dim_lambda] = sample_mutate_rand(t->curr_rand[s_dim_lambda],  points_rand(rt.points, tid), 0.1f); // ~30nm
    for(int k=0;k<curr_path->length;k++)
    {
      const int beg = curr_path->v[k].rand_beg;
      if((curr_path->v[k].mode & s_emit) && (k==0))
      {
        // k==0 makes sure next event estimation doesn't trigger this. bdpt always hands us the light path first.
        // leave s_dim_envmapvsarea and s_dim_lightsource
        t->tent_rand[beg + s_dim_envmapvsarea] = t->curr_rand[beg + s_dim_envmapvsarea];
        t->tent_rand[beg + s_dim_lightsource]  = t->curr_rand[beg + s_dim_lightsource];

        t->tent_rand[beg + s_dim_light_x] = sample_mutate_rand(t->curr_rand[beg + s_dim_light_x], points_rand(rt.points, tid), .005f);
        t->tent_rand[beg + s_dim_light_y] = sample_mutate_rand(t->curr_rand[beg + s_dim_light_y], points_rand(rt.points, tid), .005f);
        t->tent_rand[beg + s_dim_edf_x]   = sample_mutate_rand(t->curr_rand[beg + s_dim_edf_x], points_rand(rt.points, tid), .001f);
        t->tent_rand[beg + s_dim_edf_y]   = sample_mutate_rand(t->curr_rand[beg + s_dim_edf_y], points_rand(rt.points, tid), .001f);
      }
      else if(curr_path->v[k].mode & s_sensor)
      {
        t->tent_rand[beg + s_dim_image_x]    = sample_mutate_rand(t->curr_rand[beg + s_dim_image_x], points_rand(rt.points, tid), 0.005f);
        t->tent_rand[beg + s_dim_image_y]    = sample_mutate_rand(t->curr_rand[beg + s_dim_image_y], points_rand(rt.points, tid), 0.005f*view_width()/(float)view_height());
        t->tent_rand[beg + s_dim_aperture_x] = sample_mutate_rand(t->curr_rand[beg + s_dim_aperture_x], points_rand(rt.points, tid), 0.3f);
        t->tent_rand[beg + s_dim_aperture_y] = sample_mutate_rand(t->curr_rand[beg + s_dim_aperture_y], points_rand(rt.points, tid), 0.3f);
      }
      else
      {
        // leave s_dim_scatter_mode and s_dim_russian_r
        t->tent_rand[beg + s_dim_scatter_mode] = t->curr_rand[beg + s_dim_scatter_mode];
        // t->tent_rand[curr_path->v[k].rand_beg + s_dim_scatter_mode] = sample_mutate_rand(t->curr_rand[curr_path->v[k].rand_beg + s_dim_scatter_mode], points_rand(rt.points, tid), 0.01f);
        t->tent_rand[beg + s_dim_russian_r]    = t->curr_rand[beg + s_dim_russian_r];

        t->tent_rand[beg + s_dim_free_path] = sample_mutate_rand(t->curr_rand[beg + s_dim_free_path], points_rand(rt.points, tid), 0.001f);
        // do the rest uniformly, includes omega and maybe next event estimation in a collapsed (path_pop) vertex:
        for(int i=beg + s_dim_omega_x;i<beg + curr_path->v[k].rand_cnt;i++)
        {
          if(curr_path->v[k].mode & s_volume)
            t->tent_rand[i] = sample_mutate_rand(t->curr_rand[i], points_rand(rt.points, tid), 0.001f);
          else
            t->tent_rand[i] = sample_mutate_rand(t->curr_rand[i], points_rand(rt.points, tid), 0.01f);
        }
      }
    }
#endif
  }
  t->tent_throughput = 0.0f;
  t->tent_throughput_mapped = 0.0f;
  t->tent_num_contribs = 0;
  path_init(tent_path, tent_path->index, tent_path->sensor.camid);
  // only respect pixel request for large steps
  if(t->large_step && x >= 0)
    path_set_pixel(tent_path, x + pointsampler(tent_path, s_dim_image_x), y + pointsampler(tent_path, s_dim_image_y));
  // tent_path->temperature = t->temperature;
  t->mean_temperature += t->temperature;
  sampler_create_path(tent_path);
}

void pointsampler_prepare_frame(pointsampler_t *s) {}
void pointsampler_stop_learning(pointsampler_t *s) {}

#undef POINTSAMPLER_INIT_SAMPLES

