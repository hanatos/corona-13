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

#include <string.h>

// number of spectral bins
#define RENDER_BINS 32

typedef struct render_t
{
  int overlays;
  float *pixel;   // colorspace buf for screenshot/display
  float *fb;      // spectral bins
}
render_t;

static inline render_t *render_init()
{
  render_t *r = (render_t *)malloc(sizeof(render_t));
  r->fb = (float *)malloc(sizeof(float)*rt.width*rt.height*RENDER_BINS);
  r->pixel = common_alloc(16, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(r->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*RENDER_BINS);
  memset(r->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  r->overlays = 0;
  return r;
}

static inline void render_cleanup(render_t *r)
{
  free(r->fb);
  free(r->pixel);
  free(r);
}

static inline void render_accum(const float i, const float j, const float lambda, const float p)
{
  float e[RENDER_BINS];
  if(!isfinite(p)) return;
  int bin = (int)(RENDER_BINS*(lambda - 400.0f)/300.0f);
  bzero(e, 4*RENDER_BINS);
  e[bin] = p;
  filter_accum(i, j, e, rt.render->fb, RENDER_BINS, 0, RENDER_BINS);
}

#include "sampler.h"

static inline void render_clear()
{
  memset(rt.render->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(rt.render->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*RENDER_BINS);
  rt.render->overlays = 0;
  sampler_clear();
  pointsampler_clear();
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
  const __m128 inv_o = _mm_set1_ps(rt.cam->iso/100.0f * 1.0f/rt.render->overlays);
#pragma omp parallel for default(none) shared(rt) schedule(static)
  for(int k=0;k<rt.width*rt.height;k++)
  {
    ((__m128*)(rt.render->pixel))[k] = _mm_setzero_ps();
    colorspace_cam_to_rgb(rt.render->fb + RENDER_BINS*k, rt.render->pixel + 4*k + 1);
    ((__m128*)(rt.render->pixel))[k] = _mm_mul_ps(((__m128*)(rt.render->pixel))[k], inv_o);
  }
  display_update(rt.display, rt.render->pixel);
}

static inline void render_print_info(FILE *fd)
{
  fprintf(fd, "render   : global illumination, %d spectral bins\n", RENDER_BINS);
  fprintf(fd, "           samples per pixel: %d\n", rt.render->overlays);

  const float inv_px = 1.0f/(rt.cam->width*rt.cam->height);
  float sumr = 0, sumg = 0, sumb = 0;
#pragma omp parallel for default(none) shared(rt) reduction(+:sumr) reduction(+:sumg) reduction(+:sumb) schedule(static)
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
  char tempname2[BUFSIZ];
  const __m128 inv_o = _mm_set1_ps(rt.cam->iso/100.0f * 1.0f/rt.render->overlays);
  for(int bin=0;bin<RENDER_BINS;bin++)
  {
#pragma omp parallel for default(none) shared(rt, bin) schedule(static)
    for(int k=0;k<rt.width*rt.height;k++)
    {
      ((__m128*)(rt.render->pixel))[k] = _mm_setzero_ps();
      for(int i=0;i<3;i++) rt.render->pixel[4*k + i + 1] = rt.render->fb[bin + RENDER_BINS*k];
      ((__m128*)(rt.render->pixel))[k] = _mm_mul_ps(((__m128*)(rt.render->pixel))[k], inv_o);
    }
    sprintf(tempname2, "%s_%d", tempname, bin);
    screenshot_write(tempname2, rt.render->pixel, 4);
  }
}
#endif
