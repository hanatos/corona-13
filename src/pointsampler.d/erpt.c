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

#include "corona_common.h"
#include "pathspace.h"
#include "pathspace/multichain.h"
#include "pathspace/nee.h"
#include "pathspace/tech.h"
#include "points.h"
#include "pointsampler.h"
#include "sampler.h"
#include "sampler_common.h"
#include "shader.h"
#include "spectrum.h"
#include "threads.h"
#include "view.h"
#include "camera.h"
#include "ext/halton/halton.h"

#include <stdio.h>
#include <float.h>

typedef struct pointsampler_thr_t
{
  float fake_rand[40];
  int fake_rand_enabled;
  double contribution;
  int perturb_end;
  double stats_avg_contrib;
  double stats_avg_contrib_cnt;
}
pointsampler_thr_t;

typedef struct pointsampler_t
{
  uint64_t reinit;
  pointsampler_thr_t *t;
  halton_t h;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: energy redistribution path tracer\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)calloc(1, sizeof(pointsampler_t));
  s->reinit = 0;
  s->t = calloc(rt.num_threads, sizeof(*s->t));
  // init halton points
  halton_init_random(&s->h, frame);
  return s;
}

void pointsampler_finalize(pointsampler_t *s) {}

void pointsampler_clear() { }

void pointsampler_cleanup(pointsampler_t *s)
{
  free(s);
}

float pointsampler(path_t *p, const int i)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->t[tid].fake_rand_enabled)
    return rt.pointsampler->t[tid].fake_rand[i];
  int v = p->length;
  const int end = p->v[v].rand_beg;
  // halton points
  const int dim = end + i;
  if(dim >= halton_get_num_dimensions())
  {
    // degenerate to pure random mersenne twister
    return points_rand(rt.points, tid);
  }
  else
    return halton_sample(&rt.pointsampler->h, dim, p->index);
}

#if 0 // extended lens perturbation (multi chain)
static double perturb(
    path_t *curr,
    path_t *tent,
    const int end) // last moving vertex
{
  const double a = multichain_perturb_connect(curr, tent, end);
  if(!(a > 0.0))
  {
    tent->length = 0;
    return 0.0;
  }
  return a;
}
#else // single step lens, then brownian bridge style
float smoothstep(float edge0, float edge1, float x)
{
  // Scale, bias and saturate x to 0..1 range
  x = CLAMP((x - edge0)/(edge1 - edge0), 0.0, 1.0); 
  // Evaluate polynomial
  return x*x*(3 - 2*x);
}

static double perturb(
    path_t *curr,
    path_t *tent,
    const int end) // last moving vertex
{
  *tent = *curr; // XXX remove (need to copy over vertices selectively)
  double T = 1.0;
  tent->sensor = curr->sensor;
  tent->index = curr->index;
  const int tid = common_get_threadid();
  float g1, g2;
  const float r1 = points_rand(rt.points, tid);
  const float r2 = points_rand(rt.points, tid);
  sample_gaussian(r1, r2, &g1, &g2);
  const float px = 3.f; // this is one sigma width of the jump
  tent->sensor.pixel_i += g1 * px;
  tent->sensor.pixel_j += g2 * px;
  // mutate point on aperture, after mutating outgoing direction (pixel)
  tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_set = 1;
  tent->sensor.pixel_set = 1;
  tent->lambda = curr->lambda;
  // XXX would need to eval the rest of the path, too!
  // XXX probably pays off more to do some hwl in post for the NEE portions
  // tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
  tent->time = curr->time; // <= TODO

      shader_exterior_medium(tent);
      // sensor modified from outside
      if(!(tent->sensor.pixel_i >= 0 && tent->sensor.pixel_j >= 0 &&
           tent->sensor.pixel_i < view_width() && tent->sensor.pixel_j < view_height()))
        return 0.0f;

      // sample new outgoing direction and corresponding dwp pdf.
      // G required to convert p(w_1) to p(x_1) cancels with measurement.
      // we're actually interested in the ratio of perturbation pdf, but
      // these are composed of p() * J and p() is symmetric.
      float throughput = view_cam_sample(tent);
      if(throughput <= 0.0) return 0.0f;
      T *= view_cam_pdf(curr, 0) / tent->v[1].pdf;

  T *= multichain_perturb_distance(curr, tent, 0);
  if(!(T > 0)) return 0.0;

  // now got an inited tent->v[1]

  // loop over v=2..end and lerp in the curr pos with t=(v-1)/(end+1-1)
  float delta[3] = {
    tent->v[1].hit.x[0]-curr->v[1].hit.x[0],
    tent->v[1].hit.x[1]-curr->v[1].hit.x[1],
    tent->v[1].hit.x[2]-curr->v[1].hit.x[2]};
  float dist = 0.0, dist2 = 0.0;
  for(int v=2;v<=end+1;v++)
    dist += curr->e[v].dist;
  for(int v=2;v<=end+1;v++)
  {
    float x[3];
    dist2 += curr->e[v].dist;
    // const float t = smoothstep(0.0, 1.0, 1.0-(v-1.0)/end);
    const float t = smoothstep(0.0, 1.0, 1.0-dist2/dist);
    for(int k=0;k<3;k++)
      x[k] = tent->v[v].hit.x[k] = curr->v[v].hit.x[k] + t*delta[k];

    // test connection
    if(path_project(tent, v, s_propagate_mutate) ||
        (tent->v[v].flags           != curr->v[v].flags) ||
        (tent->v[v].hit.shader      != curr->v[v].hit.shader) ||
        (tent->v[v].interior.shader != curr->v[v].interior.shader) ||
        (primid_invalid(tent->v[v].hit.prim) != primid_invalid(curr->v[v].hit.prim))) 
      return 0.0f;
    // check whether we actually arrived at vertex c
    for(int k=0;k<3;k++)
      if(fabsf(tent->v[v].hit.x[k] - x[k]) > HALFVEC_REL_SPATIAL_EPS *
          MAX(MAX(fabsf(tent->v[v].hit.x[k]), fabsf(x[k])), 1.0))
        return 0.0f;

    // TODO: double check the acceptance probability
    // if(v == 2)
    // {
    T *= shader_brdf(tent, v-1);// * shader_vol_transmittance(tent, v) * path_G(tent, v);
    T /= shader_brdf(curr, v-1);// * shader_vol_transmittance(curr, v) * path_G(curr, v);
    // }
    T *= shader_vol_transmittance(tent, v) * path_G(tent, v);
    T /= shader_vol_transmittance(curr, v) * path_G(curr, v);
    // fprintf(stderr, "bsdf: %g %g\n", shader_brdf(tent, v-1), shader_brdf(curr, v-1));
    if(!(T >= 0.0)) return 0.0;
  }
  if(end+1 == tent->length-1)
  {
    T *= lights_eval_vertex(tent, tent->length-1);
    T /= lights_eval_vertex(curr, curr->length-1);
  }
  else
  {
    T *= shader_brdf(tent, end+1);
    T /= shader_brdf(curr, end+1);
  }
  // nan-check:
  if(!(T >= 0.0)) return 0.0;
  return T;
}
#endif

static void explore(path_t *path, float value)
{
#if 0 // XXX
  view_splat(path, value);
  return;
#endif // XXX
  const int chains = 10;
  const int mutations = 1;
  path_t pdata0, pdata1;

  // TODO: multiple-try mlt says the acceptance probability should be
  // sum w(curr|path_j) / sum w(tent|path_j)
  // and there would be another set of competitor samples and complicated.
  const int tid = common_get_threadid();

  FILE *f = 0;
  int32_t cnt = 0;
  if(0)//path->length > 100 && tid == 0)
  {
    f = fopen("/tmp/path.obj", "w");
    fprintf(f, "o base_path\n");
    for(int k=0;k<path->length;k++)
      fprintf(f, "v %g %g %g\n",
          path->v[k].hit.x[0],
          path->v[k].hit.x[1],
          path->v[k].hit.x[2]);
    for(int k=1;k<path->length;k++)
      fprintf(f, "f %d %d\n", k, k+1);
    cnt = path->length;
  }

  const int end = rt.pointsampler->t[tid].perturb_end;
  double mean_a = 0.0;
  double w_center = 0.0;
  for(int c=0;c<chains;c++)
  {
    path_t *tent = &pdata0;
    path_t *curr = path;
    double w_tent = 0.0;
    double w_curr = 0.0;
    for(int m=0;m<mutations;m++)
    {
      double T = perturb(curr, tent, end);

      const double a = 1;// XXX CLAMP(T, 0.0, 1.0);
      if(!(a > 0)) continue;
      w_tent = a * value / (chains*mutations);
      w_curr += (1.0-a) * value / (chains*mutations);

      mean_a += a;
      if(f)
      {
        fprintf(f, "o chain_%03d_mut_%03d\n", c, m);
        for(int k=0;k<=end+1;k++)
          fprintf(f, "v %g %g %g\n",
              tent->v[k].hit.x[0],
              tent->v[k].hit.x[1],
              tent->v[k].hit.x[2]);
        for(int k=1;k<=end+1;k++)
          fprintf(f, "f %d %d\n", cnt+k, cnt+k+1);
        cnt += end+2;
      }
      if(points_rand(rt.points, tid) < a)
      { // accept and swap
        render_splat(curr, w_curr);
        w_curr = w_tent;
        if(curr != path)
        {
          path_t *tmp = curr;
          curr = tent;
          tent = tmp;
        }
        else
        {
          curr = tent;
          tent = &pdata1;
        }
      }
      else
      {
        render_splat(tent, w_tent);
      }
    }
    // accumulate rest of current state of markov chain:
    render_splat(curr, w_curr);
    // if(mutations > 1)
    // else w_center += w_curr;
  }
  // accumulate the center path once in case we only do multi-chains and no multi-mutations:
  // XXX if(mutations == 1) render_splat(path, w_center);
  // fprintf(stderr, "mean a %g center contrib %g\n", mean_a/(mutations*chains), w_center);
  rt.pointsampler->t[tid].stats_avg_contrib += mean_a/(mutations*chains);
  rt.pointsampler->t[tid].stats_avg_contrib_cnt ++;
  // splat accumulated rejected states:
  // in particular that means we should want high acceptance, 1-a < 1/chains
  // TODO: write out values of acceptance and center_contrib, tweak to make them about equal?
  // render_splat(path, center_contrib*value/chains);

  if(f)
  {
    fclose(f);
    // bam!
    exit(0);
  }
}

void pointsampler_splat(path_t *path, float value)
{
  if(!(value > 0.0 && value < FLT_MAX)) return;

  const int end = rt.pointsampler->t[common_get_threadid()].perturb_end;
  if(path->length <= end+1)
  {
    render_splat(path, value);
    return;
  }

  const int vk = path->length-1;
  if(path->v[vk].tech == s_tech_nee && path->length == end+2)
  {
    render_splat(path, value);
    return;
  }

  // store contribution for redistribution:
  rt.pointsampler->t[common_get_threadid()].contribution += value;
}

int pointsampler_accept(path_t *curr, path_t *tent)
{
#if 0 // i think c > 0 below is enough to be tested.
  const int end = rt.pointsampler->t[common_get_threadid()].perturb_end;
  if(tent->length <= end+1) return 0;
  const int vk = tent->length-1;
  if(tent->v[vk].tech == s_tech_nee && tent->length == end+2) return 0;
#endif

  // explore these paths:
  // collected values for splatting in poinsampler_splat() for common prefix.
  // so now we can explore one example path up to this point and pretend erpt
  // had done the exact same thing for all these paths. this is just an
  // optimisation, we do the exploration once and redistribute the sum of
  // contributions directly (instead of the expensive naive sum with multiple
  // explorations in pointsampler_splat() directly)
  double c = rt.pointsampler->t[common_get_threadid()].contribution;
  if(c > 0) explore(tent, c);
  return 0;
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  rt.pointsampler->t[common_get_threadid()].contribution = 0;
  rt.pointsampler->t[common_get_threadid()].perturb_end = 2;
  path_init(tent, tent->index, tent->sensor.camid);
  sampler_create_path(tent);
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  rt.pointsampler->t[common_get_threadid()].contribution = 0;
  rt.pointsampler->t[common_get_threadid()].perturb_end = 2;
  path_init(tent, tent->index, tent->sensor.camid);
  path_set_pixel(tent, i, j);
  sampler_create_path(tent);
}

void pointsampler_enable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->t[tid].fake_rand_enabled = 1;
}

void pointsampler_disable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->t[tid].fake_rand_enabled = 0;
}

void pointsampler_set_fake_random(pointsampler_t *s, int dim, float rand)
{
  const int tid = common_get_threadid();
  s->t[tid].fake_rand[dim] = rand;
}

void pointsampler_prepare_frame(pointsampler_t *s)
{
  double avgc = 0.0;
  for(int k=0;k<rt.num_threads;k++)
  {
    avgc += s->t[k].stats_avg_contrib/(rt.num_threads*s->t[k].stats_avg_contrib_cnt);
  }
  fprintf(stderr, "\r[pointsampler] avg contrib %g  ", avgc);
  // would wrap around int limit and run out of bits for radical inverse?  stop
  // to re-init the random bit permutations (we only pass a seed, will be
  // randomised internally):
  if(rt.threads->end >= (1ul<<32))
    halton_init_random(&s->h, rt.anim_frame + ++s->reinit);
}
