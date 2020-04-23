/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RENDER_H
#define RENDER_H

#include "corona_common.h"
#include "camera.h"
#include "display.h"
#include "points.h"
#include "spectrum.h"
#include "sampler_common.h"
#define FILTER_ATOMIC   // needed for spline filter even in the absence of lt :(
#include "filter.h"
#include "screenshot.h"

#include "sampler_common.h"

#include <string.h>

/** speed of light in meters per second. in corona-6, 10.0f = 1m in world space. */
static const double render_speed_of_light = 299792458.;
/** 33. pico seconds for one centimeter (= 0.1f resolution. */
static const double render_time_for_1cm   = .00000000003335640951; // 33.35640951e-12;
/** how far does light travel in 1 pico second. */
// static const double render_dist_m_1ps       = 10e-12 * render_speed_of_light;

typedef struct render_t
{
  int overlays;
  float *pixel;   // colorspace buf for screenshot/display
  float **fb;     // time bins

  int   time_bins;
  float time_resolution;
  float time_start;

  // laser positions, directions, wavelengths, and spectral radiances
  int ls_cnt;
  float *ls_pos;
  float *ls_dir;
  float *ls_lambda;
  float *ls_L;
}
render_t;

static inline render_t *render_init()
{
  render_t *r = (render_t *)malloc(sizeof(render_t));

  r->time_bins = 100;
  r->time_resolution = 10. * 10e-12 * render_speed_of_light; // distance in dm of 1 pico sec
  r->time_start = 0.0f; // start at t = 0 => shortest direct path from laser to camera
  r->ls_cnt = 0;

  // parse config file:
  FILE *f = fopen("temporal.cfg", "rb");
  if(f)
  {
    int rv = 0;
    rv = fscanf(f, "bins %d\n", &r->time_bins);
    assert(rv == 1);
    rv = fscanf(f, "res %f\n", &r->time_resolution);
    assert(rv == 1);
    rv = fscanf(f, "start %f\n", &r->time_start);
    assert(rv == 1);
    rv = fscanf(f, "lasers %d\n", &r->ls_cnt);
    assert(rv == 1);
    r->ls_pos = (float *)malloc(3*sizeof(float)*r->ls_cnt);
    r->ls_dir = (float *)malloc(3*sizeof(float)*r->ls_cnt);
    r->ls_lambda = (float *)malloc(sizeof(float)*r->ls_cnt);
    r->ls_L = (float *)malloc(sizeof(float)*r->ls_cnt);
    for(int k=0;k<r->ls_cnt;k++)
    {
      rv = fscanf(f, "%f %f %f %f %f %f %f %f\n", r->ls_pos+3*k+0, r->ls_pos+3*k+1, r->ls_pos+3*k+2,
                                                 r->ls_dir+3*k+0, r->ls_dir+3*k+1, r->ls_dir+3*k+2,
                                                 r->ls_lambda+k, r->ls_L+k);
      assert(rv == 8);
      normalise(r->ls_dir + 3*k);
    }
    printf("[render_temporal_init] loaded %d laser light sources\n", r->ls_cnt);
    printf("[render_temporal_init] rendering %d bins with %f dm spatial resolution\n", r->time_bins, r->time_resolution);
    printf("[render_temporal_init] starting recording after the light travelled %f dm\n", r->time_start);
    fclose(f);
  }
  else
  {
    fprintf(stderr, "[render_temporal_init] could not parse file `temporal.cfg'!\n");
  }

  r->fb = (float **)malloc(sizeof(float**)*r->time_bins);
  for(int k=0;k<r->time_bins;k++)
  {
    r->fb[k] = (float *)malloc(sizeof(float)*rt.width*rt.height*3);
    memset(r->fb[k], 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  }
  r->pixel = common_alloc(16, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(r->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  r->overlays = 0;
  return r;
}

static inline void render_cleanup(render_t *r)
{
  for(int k=0;k<r->time_bins;k++) free(r->fb[k]);
  free(r->fb);
  free(r->ls_pos);
  free(r->ls_dir);
  free(r->ls_lambda);
  free(r->ls_L);
  free(r->pixel);
  free(r);
}

static inline void render_accum(const float i, const float j, const float lambda, const float p) {}

#include "sampler.h"
#ifndef SAMPLER_DUMMY_H
  #error "render temporal requires the dummy sampler module!"
#endif
#define SAMPLER_NUM_DIRECT_LENS 1
// we're ignoring most sampler functionality (light starts at camera)
static inline sampler_t *sampler_init() {return NULL;}
static inline void sampler_cleanup(sampler_t *s) {}
static inline void sampler_clear() {}
static inline void sampler_init_light(const unsigned int prim, const unsigned int num_prims, const float L, const float k) {}
static inline void sampler_prepare_frame() {}


static inline void render_accum_dist(const float i, const float j, const float lambda, const float p, const float dist)
{
  render_t *r = rt.render;
  const int k = (int)fminf(r->time_bins-1, (dist-r->time_start) / r->time_resolution);
  float e[3];
  if(!isfinite(p)) return;
  // spectrum_p_to_cam(lambda, fminf(p, 10.0f), e);
  spectrum_p_to_cam(lambda, p, e);
  filter_accum(i, j, e, rt.render->fb[k], 3, 0, 3);
}

static inline void render_clear()
{
  memset(rt.render->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  for(int k=0;k<rt.render->time_bins;k++)
    memset(rt.render->fb[k], 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  rt.render->overlays = 0;
  sampler_clear();
  pointsampler_clear();
}

static inline void sampler_create_paths(const float unused_i, const float unused_j)
{
  ray_t ray;
  rayhit_t hit, hitl;
  pointsampler_thr_t *t = rt.pointsampler->t + common_get_threadid();
  hit.adjoint = hitl.adjoint = 0;
  hit.prim = hitl.prim = -1;
  float contrib, sensor, dist = 0.0f, dist_offs = 0.0f;
  float brdf, omega_in[3];
  float i = 0, j = 0;

  const int ls = t->rand[PS_LAMBDA] * rt.render->ls_cnt;
  contrib = rt.render->ls_L[ls];
  for(int k=0;k<3;k++) ray.pos[k] = rt.render->ls_pos[k];
  for(int k=0;k<3;k++) ray.dir[k] = rt.render->ls_dir[k];
  for(int k=0;k<3;k++) ray.invdir[k] = 1.0/rt.render->ls_dir[k];
  rayhit_init(&hit);
  hitl.lambda = hit.lambda = rt.render->ls_lambda[ls];
  hit.dist = INFINITY;
  hit.adjoint = 0;

  while(1)
  {
    accel_intersect(rt.accel, &ray, &hit);

    if(hit.prim == -1) return; // lost the sample

    if(dist == 0.0f) dist_offs = hit.dist;
    dist += hit.dist;
    // if((int)(dist / RENDER_RESOLUTION) >= RENDER_NUM_TIME_BINS) return;
    for(int k=0;k<3;k++) omega_in[k] = ray.dir[k];
    float rr = POINTSAMPLER_RAND(1);
    float rx = POINTSAMPLER_RAND(2);
    float ry = POINTSAMPLER_RAND(3);
    float ra = points_rand(rt.points, common_get_threadid());
    contrib *= shader_prepare(&ray, &hit, POINTSAMPLER_RAND(0), rr);

    // deterministically connect to lens (hit prim might be lost in medium)
    if(hit.prim >= 0) for(int l=0;l<SAMPLER_NUM_DIRECT_LENS;l++)
    {
      // get hitl on lens randomly (camera_sample_lens), clip to viewing frustum
      float x   = POINTSAMPLER_RANDN(PS_LENS_X),  // points_rand(rt.points, threadid),
            y   = POINTSAMPLER_RANDN(PS_LENS_Y);  // points_rand(rt.points, threadid);
      hitl.lambda = hit.lambda;
      sensor = camera_get_pixel(rt.cam, &hit, &i, &j, &hitl, x, y);
      if(sensor == 0.0f) continue;
      // if(prim_get_ray(&hit, &hitl, &ray)) continue;
      // if(dotproduct(hitl.normal, ray.dir) >= 0) continue;
      if(prim_get_ray(&hitl, &hit, &ray)) continue;
      if(dotproduct(hitl.normal, ray.dir) <= 0) continue;
      // if(dotproduct(hit.normal, ray.dir) <= 0) continue; // inside glass, cloth ..?
      if(accel_visible(rt.accel, &ray))
      {
        const float r2 = fminf(FLT_MAX/2.0f, 1.0f/dotproduct(ray.dir, ray.dir));
        // for(int k=0;k<3;k++) ray.dir[k] *= sqrtf(r2);
        for(int k=0;k<3;k++) ray.dir[k] *= - sqrtf(r2);
        float brdf = shader_brdf(omega_in, &hit, ray.dir);
        //      const float g = fmaxf(0.0f, - dotproduct(hitl.normal, ray.dir))*fmaxf(0.0f, dotproduct(hit.normal, ray.dir))*r2;
        //const float g = - dotproduct(hitl.normal, ray.dir)*fabsf(dotproduct(hit.normal, ray.dir))*r2;
        const float g = - dotproduct(hitl.normal, ray.dir)*fabsf(dotproduct(hit.normal, ray.dir))*r2;
        const float adj = fmaxf(0.0f, dotproduct(hit.normal, omega_in)/dotproduct(hit.normal, omega_in));
        brdf *= sensor*contrib * g * adj / (float)SAMPLER_NUM_DIRECT_LENS;
        // pointsampler_accum(i, j, hit.lambda, brdf);
        render_accum_dist(i, j, hit.lambda, brdf, dist + 1.0/sqrtf(r2) - dist_offs);
      }
    }
    brdf = shader_sample(omega_in, &hit, &ray, rx, ry, rr);

    contrib *= brdf;
    if(t->pathlen++ > 4)
    {
      // russian roulette
      const float thrs = fminf(brdf, 0.9);
      if(ra >= thrs) return;
      contrib *= 1.0f/thrs;
    }
    if(t->pathlen > POINTSAMPLER_MAX_PATH_LENGTH || brdf == 0.0f) return;
  }
}

static inline void render()
{
#pragma omp parallel for default(none) shared(rt) schedule(dynamic)
  for(int p=0;p<rt.width*rt.height;p++)
  {
    // for first few, new sample
    pointsampler_next_sample(rt.pointsampler);
    sampler_create_paths(p%rt.width, p/rt.width);
    pointsampler_add_sample(rt.pointsampler);
    // get next point (and use its components for next sample)
    points_next(rt.points, common_get_threadid());
  }
  rt.render->overlays++;
#ifndef DISPLAY_NULL_H
  const __m128 inv_o = _mm_set1_ps(rt.cam->iso/100.0f * 1.0f/rt.render->overlays);
#pragma omp parallel for default(none) shared(rt) schedule(static)
  for(int k=0;k<rt.width*rt.height;k++)
  {
    ((__m128*)(rt.render->pixel))[k] = _mm_setzero_ps();
    for(int i=0;i<rt.render->time_bins;i++)
    {
      __m128 tmp;
      colorspace_cam_to_rgb(rt.render->fb[i] + 3*k, ((float *)&tmp)+1);
      ((__m128*)(rt.render->pixel))[k] = _mm_add_ps(((__m128*)(rt.render->pixel))[k], _mm_mul_ps(tmp, inv_o));
    }
  }
#endif
  display_update(rt.display, rt.render->pixel);
}

static inline void render_print_info(FILE *fd)
{
  fprintf(fd, "render   : temporal ray tracing\n");
  fprintf(fd, "           samples per pixel: %d\n", rt.render->overlays);
  filter_print_info(fd);
}

static inline void render_screenshot(const char *tempname)
{
  char filename[1024];
  const __m128 inv_o = _mm_set1_ps(rt.cam->iso/100.0f * 1.0f/rt.render->overlays);
  for(int i=0;i<rt.render->time_bins;i++)
  {
#pragma omp parallel for default(none) shared(rt, i) schedule(static)
    for(int k=0;k<rt.width*rt.height;k++)
    {
      ((__m128*)(rt.render->pixel))[k] = _mm_setzero_ps();
      __m128 tmp;
      colorspace_cam_to_rgb(rt.render->fb[i] + 3*k, ((float *)&tmp)+1);
      ((__m128*)(rt.render->pixel))[k] = _mm_add_ps(((__m128*)(rt.render->pixel))[k], _mm_mul_ps(tmp, inv_o));
    }
    snprintf(filename, 1024, "%s_time%04d", tempname, i);
    screenshot_write(filename, rt.render->pixel, 4);
  }
}

#endif
