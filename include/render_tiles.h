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

typedef struct blocks_t
{
  int32_t size;
  int32_t offx, offy;
  int32_t x,y,dx,dy;
  int32_t t;
  int32_t i, max_i;
}
blocks_t;

static inline void render_blocks_init(blocks_t *b, const int tilewd, const int wd, const int ht)
{
  b->size = tilewd;
  b->offx = (wd-1)/tilewd/2;
  b->offy = (ht-1)/tilewd/2;
  b->t = MAX((wd/tilewd), (ht/tilewd));
  b->max_i = b->t*b->t;
}

static inline void render_blocks_reset(blocks_t *b)
{
  b->x = b->y = b->dx = 0;
  b->dy = -1;
  b->i = 0;
}

static inline void render_blocks_step(blocks_t *b)
{
  if((b->x == b->y) ||
      ((b->x < 0) && (b->x == -b->y)) ||
      ((b->x > 0) && (b->x == 1-b->y)))
  {
    b->t = b->dx;
    b->dx = -b->dy;
    b->dy = b->t;
  }
  b->x += b->dx;
  b->y += b->dy;
  b->i++;
}

static inline int render_blocks_get_offsets(
    const int tilewd,
    const int wd,
    const int ht,
    int *x,
    int *y)
{
  blocks_t b;
  render_blocks_init(&b, tilewd, wd, ht);
  render_blocks_reset(&b);
  int i = 1;
  x[0] = tilewd*(b.x+b.offx);
  y[0] = tilewd*(b.y+b.offy);
  while(1)
  {
    render_blocks_step(&b);
    // only count tiles inside the image. assume first block is inside.
    if(b.x + b.offx >= 0 && b.y + b.offy >= 0 && tilewd*(b.x+b.offx) + tilewd <= rt.width && tilewd*(b.y+b.offy) + tilewd <= rt.height)
    {
      x[i] = tilewd*(b.x+b.offx);
      y[i] = tilewd*(b.y+b.offy);
      i++;
    }
    if(b.i == b.max_i) return i;
  }
}


typedef struct render_t
{
  int overlays;      // how many full frames rendered out so far

  int num_tiles;
  int tile;          // global tile number work queue

  int *offx, *offy;  // pixel offset of tiles
  int *tile_index;   // tile being worked on by thread id
  int *pixel_index;  // pixel being worked on by thread id

  float *pixel;      // colorspace buf for screenshot/display
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
  memset(&r->path0, 0, sizeof(path_t));
  memset(&r->path1, 0, sizeof(path_t));
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
  r->pixel = common_alloc(16, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(r->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  memset(r->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);

  const int tilewd = 32;
  r->tile = 0;
  r->num_tiles = (rt.width/tilewd)*(rt.height/tilewd);
  r->offx = (int*)malloc(sizeof(int)*r->num_tiles);
  r->offy = (int*)malloc(sizeof(int)*r->num_tiles);
  r->tile_index = (int*)malloc(sizeof(int)*rt.num_threads);
  r->pixel_index = (int*)malloc(sizeof(int)*rt.num_threads);
  int num_tiles = render_blocks_get_offsets(tilewd, rt.width, rt.height, r->offx, r->offy);
  assert(num_tiles == r->num_tiles);

  r->overlays = 0;
  r->time_wallclock = common_time_wallclock();
  r->time_user = common_time_user();
  return r;
}

static inline void render_cleanup(render_t *r)
{
  free(r->offx);
  free(r->offy);
  free(r->tile_index);
  free(r->pixel_index);
  free(r->fb);
  free(r->pixel);
  free(r);
}

static inline void render_clear()
{
  memset(rt.render->pixel, 0, sizeof(float)*rt.cam->width*rt.cam->height*4);
  memset(rt.render->fb, 0, sizeof(float)*rt.cam->width*rt.cam->height*3);
  rt.render->overlays = 0;
  rt.render->tile = 0;
  memset(rt.render->tile_index, 1, sizeof(int)*rt.num_threads);
  memset(rt.render->pixel_index, 1, sizeof(int)*rt.num_threads);
  sampler_clear(rt.sampler);
  pointsampler_clear();
  rt.render->time_wallclock = common_time_wallclock();
  rt.render->time_user = common_time_user();
}

void render_prepare_frame();
void render_accum(path_t *p, const float value);
void render_sample(uint64_t sample);
void render_update(uint64_t pixel);

static inline void render_print_info(FILE *fd)
{
  fprintf(fd, "render   : global illumination, memory friendly framebuffer\n");
  fprintf(fd, "           samples per pixel: %d\n", rt.render->overlays);
  fprintf(fd, "           elapsed times: wallclock %.1f s user %.1f s\n",
      common_time_wallclock() - rt.render->time_wallclock,
      common_time_user() - rt.render->time_user);

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
