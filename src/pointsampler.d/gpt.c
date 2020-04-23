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
#include "pointsampler.h"
#include "sampler.h"
#include "render.h"
#include "lights.h"
#include "pathspace.h"
#include "points.h"
#include "spectrum.h"
#include "shader.h"
#include "pathspace/tech.h"
#include "threads.h"
#include "pathspace/multichain.h"
#include "ext/halton/halton.h"

#include <stdio.h>
#include <float.h>
#include <immintrin.h>

#if 0 // on GCC, compile with option -mbmi2, requires Haswell or better :(
static uint64_t xy_to_morton (uint32_t x, uint32_t y)
{
  return _pdep_u32(x, 0x55555555) | _pdep_u32(y,0xaaaaaaaa);
}

static void morton_to_xy (uint64_t m, uint32_t *x, uint32_t *y)
{
  *x = _pext_u64(m, 0x5555555555555555);
  *y = _pext_u64(m, 0xaaaaaaaaaaaaaaaa);
}
#endif

#if 0
static uint32_t morton_1(uint32_t x)
{
  // x = x & 0x5555555555555555;
  // x = (x | (x >> 1)) & 0x3333333333333333;
  // x = (x | (x >> 2)) & 0x0F0F0F0F0F0F0F0F;
  // x = (x | (x >> 4)) & 0x00FF00FF00FF00FF;
  // x = (x | (x >> 8)) & 0x0000FFFF0000FFFF;
  // x = (x | (x >> 16)) & 0xFFFFFFFFFFFFFFFF;
  x = x & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0F0F0F0F;
  x = (x | (x >> 4)) & 0x00FF00FF;
  x = (x | (x >> 8)) & 0x0000FFFF;
  return x;
}

static void morton_to_xy(uint64_t m, uint32_t *x, uint32_t *y)
{
  *x = morton_1(m);
  *y = morton_1(m >> 1);
}
#endif

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
  halton_t h;
  uint64_t reinit;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: halton points\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  s->reinit = 0;
  halton_init_random(&s->h, frame);
  return s;
}

int pointsampler_accept(path_t *curr, path_t *tent) { return 0; }
void pointsampler_clear() {}
void pointsampler_cleanup(pointsampler_t *s)
{
  free(s->rand);
  free(s);
}
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}
void pointsampler_finalize(pointsampler_t *s) {}

float pointsampler(path_t *p, int i)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[i];

  int v = p->length;
  const int end = p->v[v].rand_beg;
  const int dim = end + i;
  if(dim >= halton_get_num_dimensions())// || dim == s_dim_lambda) // XXX choose s_dim_lambda random, too
    // degenerate to pure random mersenne twister
    // XXX may need to fix those to actually produce deterministic results!
    return points_rand(rt.points, common_get_threadid());
  else
    // note that this clips the bits in p->index to 32:
    // XXX moved up by two to reuse s_dim_image_{x,y} XXX also s_dim_lambda
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
    const int shift_x,
    const int shift_y,
    const int end) // last moving vertex (i.e. end==1 means vpl style reconnection)
{
  path_copy(tent, curr);
  double T = 1.0;
  // tent->sensor = curr->sensor;
  // tent->index = curr->index;
  // const int tid = common_get_threadid();
  // float g1, g2;
  // const float r1 = points_rand(rt.points, tid);
  // const float r2 = points_rand(rt.points, tid);
  // sample_gaussian(r1, r2, &g1, &g2);
  // const float px = 3.f; // this is one sigma width of the jump
  // tent->sensor.pixel_i += g1 * px;
  // tent->sensor.pixel_j += g2 * px;
  tent->sensor.pixel_i += shift_x;
  tent->sensor.pixel_j += shift_y;
  // mutate point on aperture, after mutating outgoing direction (pixel)
  // tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.3f);
  // tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_set = 1;
  tent->sensor.pixel_set = 1;
  // tent->lambda = curr->lambda;
  // XXX would need to eval the rest of the path, too!
  // XXX probably pays off more to do some hwl in post for the NEE portions
  // tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
  // tent->time = curr->time; // <= TODO

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

  if(tent->length == 2)
  {
    T *= lights_eval_vertex(tent, tent->length-1);
    T /= lights_eval_vertex(curr, curr->length-1);
    return T;
  }

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
void pointsampler_splat(path_t *p, float value)
{
  // XXX TODO multiplex into more buffers: primal, primal+x, primal+y
  // XXX TODO: see how we could do this in accept() only, not for every splat
#if 1
  render_splat(p, value); // regular mono wavelength version
  path_t mut;
  double Jp, Jm;
  double wp, wm;
  
  Jm = perturb(p, &mut, -1,  0, 1);
  Jp = perturb(p, &mut,  1,  0, 1);
  wp = Jp > 0 ? 1.0 : 0.0; // TODO MIS
  wm = Jm > 0 ? 1.0 : 0.0;
  wp /= wm+wp;
  wm /= wm+wp;
  mut.sensor = p->sensor;
  mut.sensor.camid = 1;
  if(wm > 0) render_splat(&mut, wm*(Jm*value-value));
  if(wp > 0) render_splat(&mut, wp*(Jp*value-value));

  Jm = perturb(p, &mut, 0, -1, 1);
  Jp = perturb(p, &mut, 0,  1, 1);
  wp = Jp > 0 ? 1.0 : 0.0;
  wm = Jm > 0 ? 1.0 : 0.0;
  wp /= wm+wp;
  wm /= wm+wp;
  mut.sensor = p->sensor;
  mut.sensor.camid = 2;
  if(wm > 0) render_splat(&mut, wm*(Jm*value-value));
  if(wp > 0) render_splat(&mut, wp*(Jp*value-value));
#else
  // splat with a couple hero wavelegths
#ifdef HEROWAVELENGTHS
  const int num_wavelengths = HEROWAVELENGTHS;
#else
  const int num_wavelengths = 4;
#endif
  double pdfsum = 0.0f;
  double val[num_wavelengths];
  double pdf[num_wavelengths];
  double mis[num_wavelengths];
  double lam[num_wavelengths];
  double shared_val0 = 1.0;
  path_t perpath;

  for(int i=0;i<num_wavelengths;i++)
  {
    path_copy(&perpath, p);
#if 1
    float l_rand = (p->lambda - spectrum_sample_min)/(spectrum_sample_max - spectrum_sample_min);
    l_rand = fmodf(l_rand + i/(float)num_wavelengths, 1.0f);
#else
    float l_rand = points_rand(rt.points, common_get_threadid());
#endif
    perpath.lambda = lam[i] = spectrum_sample_min + (spectrum_sample_max - spectrum_sample_min)*l_rand;

    for(int v=1; v<p->length; v++)
    {
      if(path_edge_init_volume(&perpath, v)) goto fail;
      // catch volume stack propagation edge case:
      // TODO: make sure these paths aren't passed on until here (inconsistent mode coming from measurement_contrib eval vs sample)
      if(!(perpath.v[v].flags & s_environment) && !(perpath.v[v].mode & s_sensor) &&
          primid_invalid(perpath.v[v].hit.prim) && (perpath.e[v].vol.shader < 0)) goto fail;
      if(!((perpath.v[v].mode & s_volume) && shader_vol_hete(&perpath, v))) // XXX assume hete volumes are monochromatic
        shader_prepare(&perpath, v);
    }
#if 0
    val[i] = path_measurement_contribution_dx(&perpath, 0, perpath.length-1);
#else // XXX assume monochromatic media and optimise away transmittance call:
    if(i == 0)
    {
    double m = 1.0;
    for(int v=1;v<perpath.length-1;v++)
    {
      m *= perpath.e[v].transmittance;
      double brdf = shader_brdf(&perpath, v);
      m *= brdf;
      if(!((perpath.v[v].mode & s_volume) && shader_vol_hete(&perpath, v))) // XXX assume hete volumes are monochromatic
        shared_val0 *= brdf;
      m *= path_G(&perpath, v);
    }
    m *= path_G(&perpath, perpath.length-1);
    m *= perpath.e[perpath.length-1].transmittance;
    double tmp = lights_eval_vertex(&perpath, perpath.length-1);
    m *= tmp;
    shared_val0 *= tmp;
    tmp = view_cam_eval(&perpath);
    m *= tmp;
    shared_val0 *= tmp;
    val[i] = m;
    }
    else
    {
      // just recompute what we know may have changed:
      val[i] = val[0] / shared_val0;
      val[i] *= lights_eval_vertex(&perpath, perpath.length-1);
      val[i] *= view_cam_eval(&perpath);
      for(int v=1;v<perpath.length-1;v++)
        if(!((perpath.v[v].mode & s_volume) && shader_vol_hete(&perpath, v))) // XXX assume hete volumes are monochromatic
          val[i] *= shader_brdf(&perpath, v);
    }
#endif
    pdf[i] = path_tech_pdf_as_sampled(&perpath);
    mis[i] = sampler_mis_weight(&perpath);
    pdfsum += pdf[i];
    if(0)
    {
fail:
      pdf[i] = 0.0f;
      val[i] = 0.0f;
    }
  }

  float acc[3] = {0.0f};
  for(int i=0;i<num_wavelengths;i++)
  {
    perpath.lambda = lam[i];

    if(!(val[i] > 0.0) || !(pdf[i] > 0.0)) continue;
    double v = val[i] / pdfsum  * mis[i];
    float col[3];
    view_deferred_splat(&perpath, v, col);
    for(int k=0;k<3;k++) acc[k] += col[k];
  }
  view_splat_col(&perpath, acc);
#endif
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  path_init(tent, tent->index, tent->sensor.camid);
#if 0
  const uint64_t ts = 32;
  const uint64_t num_px = view_width()*view_height();
  const uint64_t frame = tent->index/num_px;
  const uint64_t px = tent->index - num_px*frame;
  const uint64_t tile_idx = px / (ts*ts);
  const uint64_t ty = tile_idx/(view_width()/ts), tx = tile_idx - ((int)ty)*(view_width()/ts);
  const uint64_t ti = px - tile_idx * ts * ts;

  uint32_t seed = tent->index * 1337;
  seed *= 0x15A4E35;
  const float r0 = (seed >> 16) / (0xffff + 1.0f);
  seed *= 0x15A4E35;
  const float r1 = (seed >> 16) / (0xffff + 1.0f);
  seed *= 0x15A4E35;
  const float r2 = (seed >> 16) / (0xffff + 1.0f);
  tent->lambda = spectrum_sample_lambda(r2, 0);

  uint32_t mx, my;
  morton_to_xy (ti, &mx, &my);
  float x = ts*tx + mx + r0;//points_rand(rt.points, common_get_threadid());
  float y = ts*ty + my + r1;//points_rand(rt.points, common_get_threadid());
#endif
  float x = view_width() * halton_sample(&rt.pointsampler->h, s_dim_image_x, tent->index);
  float y = view_height() * halton_sample(&rt.pointsampler->h, s_dim_image_y, tent->index);
// #define GRAD_Y
#ifdef GRAD_X
  x = fmodf(x+1, view_width());
#endif
#ifdef GRAD_Y
  y = fmodf(y+1, view_height());
#endif
  path_set_pixel(tent, x, y);
  sampler_create_path(tent);
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  path_init(tent, tent->index, tent->sensor.camid);
  path_set_pixel(tent, i, j);
  sampler_create_path(tent);
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

void pointsampler_prepare_frame(pointsampler_t *s)
{
  // would wrap around int limit and run out of bits for radical inverse?  stop
  // to re-init the random bit permutations (we only pass a seed, will be
  // randomised internally):
  if(rt.threads->end >> 32 > s->reinit)
    halton_init_random(&s->h, rt.anim_frame + ++s->reinit);
}
