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

#ifndef RENDER_GLOBAL_ILLUMINATION_H
#define RENDER_GLOBAL_ILLUMINATION_H

#include "corona_common.h"
#include "camera.h"
#include "display.h"
#include "points.h"
#include "sampler.h"
#include "spectrum.h"
#include "sampler_common.h"
#define FILTER_ATOMIC   // needed for spline filter even in the absence of lt :(
#include "filter.h"
#include "screenshot.h"

#include <string.h>

#define ERPT_MUTATIONS (rt.width)

typedef struct render_t
{
  int overlays;
  float *pixel;   // colorspace buf for screenshot/display
  float *fb;

  double time_wallclock;
  double time_user;
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

static inline render_tls_t *render_tls_init()
{
  render_tls_t *r = (render_tls_t *)malloc(sizeof(render_tls_t));
  path_init(&r->path0);
  path_init(&r->path1);
  r->curr_path = &r->path0;
  r->tent_path = &r->path1;
  return r;
}

static inline void render_tls_cleanup(render_tls_t *r)
{
  free(r);
}

static inline render_t *render_init()
{
  render_t *r = (render_t *)malloc(sizeof(render_t));
  r->fb = (float *)malloc(sizeof(float)*rt.width*rt.height*3);
  r->pixel = common_alloc(16, sizeof(float)*(rt.cam->width*rt.cam->height+1)*4);
  memset(r->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  memset(r->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  r->overlays = 0;
  r->time_wallclock = common_time_wallclock();
  r->time_user = common_time_user();
  return r;
}

static inline void render_cleanup(render_t *r)
{
  free(r->fb);
  free(r->pixel);
  free(r);
}

static inline void render_clear()
{
  memset(rt.render->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(rt.render->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  rt.render->overlays = 0;
  sampler_clear(rt.sampler);
  pointsampler_clear();
  rt.render->time_wallclock = common_time_wallclock();
  rt.render->time_user = common_time_user();
}

static inline void render_prepare_frame() {}
void render_accum(const path_t *p, const float value);
void render_sample(uint64_t sample);
void render_update(uint64_t pixel);

static inline void _render_print_time(FILE *fd, double time)
{
  int seconds = (int)time;
  int milli = (time - seconds)*1000;
  int mins = seconds / 60.0;
  seconds -= mins * 60;
  int hours = mins / 60.0;
  mins -= hours * 60;
  int days = hours / 24.0;
  hours -= days * 24;
  if(days > 0) fprintf(fd, "%d days ", days);
  if(days > 0 || hours > 0) fprintf(fd, "%dh ", hours);
  if(days > 0 || hours > 0 || mins > 0) fprintf(fd, "%02d:", mins);
  fprintf(fd, "%02d.%02d", seconds, milli/10);
}

static inline void render_print_info(FILE *fd)
{
  const double wallclock = common_time_wallclock() - rt.render->time_wallclock;
  const double user = common_time_user() - rt.render->time_user;
  fprintf(fd, "render   : global illumination, energy redistribution path tracing with %d mutations\n", ERPT_MUTATIONS);
  fprintf(fd, "           samples per pixel: %d (%.2f s/prog) max path vertices %d\n", rt.render->overlays, wallclock/rt.render->overlays, PATHSPACE_MAX_VERTS);
  fprintf(fd, "           elapsed times: wallclock %.2fs (", wallclock);
  _render_print_time(fd, wallclock);
  fprintf(fd, "), user %.2fs (", user);
  _render_print_time(fd, user);
  fprintf(fd, ")\n");

  const float inv_o = rt.cam->iso/100.0f * 1.0f/rt.render->overlays;
// #pragma omp parallel for default(none) shared(rt) schedule(static)
  for(int k=0;k<rt.width*rt.height;k++)
    for(int i=0;i<3;i++) rt.render->pixel[4*k+1+i] = rt.render->fb[3*k + i] * inv_o;
  const float inv_px = 1.0f/(rt.cam->width*rt.cam->height);
  float sumr = 0, sumg = 0, sumb = 0;
// #pragma omp parallel for default(none) shared(rt) reduction(+:sumr) reduction(+:sumg) reduction(+:sumb) schedule(static)
  for(int k=1;k<4*rt.cam->width*rt.cam->height;k++)
  {
    sumr += inv_px*rt.render->pixel[k++];
    sumg += inv_px*rt.render->pixel[k++];
    sumb += inv_px*rt.render->pixel[k++];
  }
  fprintf(fd, "           average image intensity (rgb): (%f %f %f)\n", sumr, sumg, sumb);
  filter_print_info(fd);
}

static inline void render_screenshot(const char *tempname)
{
  pointsampler_finalize(rt.pointsampler);
  const __m128 inv_o = _mm_set1_ps(rt.cam->iso/100.0f * 1.0f/rt.render->overlays);
// #pragma omp parallel for default(none) shared(rt) schedule(static)
  for(int k=0;k<rt.width*rt.height;k++)
  {
    ((__m128*)(rt.render->pixel))[k] = _mm_setzero_ps();
    for(int i=0;i<3;i++) rt.render->pixel[4*k + i + 1] = rt.render->fb[3*k + i];
    ((__m128*)(rt.render->pixel))[k] = _mm_mul_ps(((__m128*)(rt.render->pixel))[k], inv_o);
  }
  screenshot_write(tempname, rt.render->pixel, 4);
}
#endif
