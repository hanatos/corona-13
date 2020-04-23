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
#include "pathspace.h"
#include "pathspace/guided.h"
#include "points.h"
#include "spectrum.h"
#include "threads.h"
#include "view.h"
#include "screenshot.h"
#include "filter.h"
#include "shader.h"
#include "ext/halton/halton.h"

#include <stdio.h>
#include <float.h>
#include <stdatomic.h>

#define HEROWAVELENGTH 0
#define DUMP_GUIDING_HIST 0
#define DUMP_PDF_PATH_INFO 0
#define PRINT_THROUPUT_HISTOGRAM 0
#define DUMP_BEAUTY 1
#define DUMP_BEAUTY_SPLITS 0
#define NUM_BEAUTY_SPLITS 7
#ifndef STOP_LEARNING
#define STOP_LEARNING 1
#endif
#ifndef CLEAR_AFTER_LEARNING
#define CLEAR_AFTER_LEARNING 1
#endif
#define LEARN_ITERATIONS 1025
#ifndef LEARN_TIME
#define LEARN_TIME 3600
#endif

static atomic_uint_fast64_t guided_samples = 0;
static atomic_uint_fast64_t guided_misses = 0;

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct variance_t
{
  uint64_t count;
  double M2;
  double mean;
}
variance_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
  halton_t h;
  uint64_t reinit;        // when we last inited new scrambling for halton points
  uint64_t halton_org;    // original halton begin of progression index
  uint64_t halton_beg;    // where uniform halton indices begin for the next frame
  uint64_t halton_end;    // where uniform indices end (first guiding index)
  float uniform_ratio;    // fraction of uniform samples to mix in
  float *guiding_sampling_hist;    // buffer to export guiding sampled histogram
  float *guiding_hist;    // buffer to export guiding path histogram
  
  float *fb_uniform_split[NUM_BEAUTY_SPLITS];
  float *fb_guiding_split[NUM_BEAUTY_SPLITS];
  float *fb_uniform;
  float *fb_guiding;
  float *fb_sequence;
  float *fb_sequence_nc;

  guided_cache_t *cache;  // guiding caches
  uint64_t num_throughput_bins;
  uint64_t *throughput_histogram;
  uint64_t *throughput_histogram_old;
  uint64_t *throughput_histogram_diff;

  int heap_size;
  int num_nb;

  render_mode_t render_mode;
  double time_stamp;
  uint64_t learn_iterations;
  uint64_t learn_paths;
  double learn_time;
  double start_time;
}
pointsampler_t;

static void splat_fb(float* fb, path_t* p, float val)
{
  // XXX api changed:
  return;
  // float col[3];
  // spectrum_p_to_camera(p->lambda, val, col);
  // const float x = p->sensor.pixel_i;
  // const float y = p->sensor.pixel_j;
  // filter_splat(x, y, col, fb, 3, 0, 3, view_width(), view_height());
}

static void write_fb(float* fb, const char* name, uint64_t norm)
{
  return; // XXX
  const float inv_o = 1.f/norm;
  char filename[1024];
  snprintf(filename, sizeof(filename), "%s%s_%s", rt.basename, rt.output_filename, name);
  screenshot_write(filename, fb, 3, 0, 3, view_width(), view_height(), inv_o);
}

void pointsampler_print_info(FILE *f)
{
  uint64_t histogram_max = 0;
  uint64_t histogram_max_diff = 0;
  const int lines = 3;
  for(uint64_t bin=0;bin<rt.pointsampler->num_throughput_bins;bin++)
  {
    if (rt.pointsampler->throughput_histogram[bin]>histogram_max)
      histogram_max = rt.pointsampler->throughput_histogram[bin];
    if (rt.pointsampler->throughput_histogram_diff[bin]>histogram_max_diff)
      histogram_max_diff = rt.pointsampler->throughput_histogram_diff[bin];
  }

  if (histogram_max <= 0)
  {
    fprintf(f, "\nhistogram max = %ld\n", histogram_max);
  }
  else
  {
    fprintf(f, "           ");
    fprintf(f, "histogram:                      ");
    fprintf(f, "           ");
    fprintf(f, "histogram diff (max %09ld): \n", histogram_max_diff);

    for(int h=lines-1;h>=0;h--)
    {
      fprintf(f, "           ");
      for(int bin=0;bin<rt.pointsampler->num_throughput_bins;bin++)
      {
        float level = rt.pointsampler->throughput_histogram[bin] / (float)histogram_max;
        float fill = level*lines-h;
        if(fill <= 0) fprintf(f, " ");
        else if(fill <= 1./8.) fprintf(f, "\u2581");
        else if(fill <= 2./8.) fprintf(f, "\u2582");
        else if(fill <= 3./8.) fprintf(f, "\u2583");
        else if(fill <= 4./8.) fprintf(f, "\u2584");
        else if(fill <= 5./8.) fprintf(f, "\u2585");
        else if(fill <= 6./8.) fprintf(f, "\u2586");
        else if(fill <= 7./8.) fprintf(f, "\u2587");
        else /*if(fill <= 8./8.)*/  fprintf(f, "\u2588");
      }
      fprintf(f, "           ");
      for(int bin=0;bin<rt.pointsampler->num_throughput_bins;bin++)
      {
        float level = rt.pointsampler->throughput_histogram_diff[bin] / (float)histogram_max_diff;
        float fill = level*lines-h;
        if (fill <= 0) fprintf(f, " ");
        else if(fill <= 1./8.) fprintf(f, "\u2581");
        else if(fill <= 2./8.) fprintf(f, "\u2582");
        else if(fill <= 3./8.) fprintf(f, "\u2583");
        else if(fill <= 4./8.) fprintf(f, "\u2584");
        else if(fill <= 5./8.) fprintf(f, "\u2585");
        else if(fill <= 6./8.) fprintf(f, "\u2586");
        else if(fill <= 7./8.) fprintf(f, "\u2587");
        else fprintf(f, "\u2588");
      }
      fprintf(f, "\n");
    }
    fprintf(f, "          1");
    fprintf(f, "|  5|    |    |  20|    |    |  ");
    fprintf(f, "          1");
    fprintf(f, "|  5|    |    |  20|    |    |  \n");
  }
  fprintf(f, "mutations: guided sampling\n");
  fprintf(f, "learning took %gs (%ld iterations, %ld guide paths)\n", 
      rt.pointsampler->learn_time, 
      rt.pointsampler->learn_iterations, 
      rt.pointsampler->learn_paths);
  guided_collect_stats(rt.pointsampler->cache, f);
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = malloc(sizeof(*s));
  s->render_mode = rm_learn;
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  s->reinit = 0;
  s->uniform_ratio = 1.0f;
  s->halton_beg = 0;
  s->halton_end = 0;
  s->halton_org = 0;
  s->num_throughput_bins = 32;
  s->throughput_histogram = calloc(s->num_throughput_bins, sizeof(uint64_t));
  s->throughput_histogram_old = calloc(s->num_throughput_bins, sizeof(uint64_t));
  s->throughput_histogram_diff = calloc(s->num_throughput_bins, sizeof(uint64_t));
  s->cache = guided_init(1<<20);
  s->guiding_sampling_hist = calloc(view_width()*view_height(), sizeof(float));
  s->guiding_hist = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_uniform = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_guiding = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_sequence = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_sequence_nc = calloc(view_width()*view_height()*3, sizeof(float));
  for (int i=0;i<NUM_BEAUTY_SPLITS;++i)
  {
    s->fb_uniform_split[i] = calloc(view_width()*view_height()*3, sizeof(float));
    s->fb_guiding_split[i] = calloc(view_width()*view_height()*3, sizeof(float));
  }
  
  s->heap_size = 4;
  s->num_nb = 10;

  halton_init_random(&s->h, frame);
  return s;
}

int pointsampler_accept(path_t *curr, path_t *tent) { return 0; }

void pointsampler_clear()
{
  rt.pointsampler->uniform_ratio = 1.0f;
  rt.pointsampler->halton_beg = 0;
  rt.pointsampler->halton_end = 0;
  rt.pointsampler->halton_org = 0;
  guided_clear(rt.pointsampler->cache);
  memset(rt.pointsampler->guiding_sampling_hist, 0, sizeof(float)*view_width()*view_height());
  memset(rt.pointsampler->guiding_hist, 0, sizeof(float)*view_width()*view_height()*3);
  memset(rt.pointsampler->fb_uniform, 0, sizeof(float)*view_width()*view_height()*3);
  memset(rt.pointsampler->fb_guiding, 0, sizeof(float)*view_width()*view_height()*3);
  memset(rt.pointsampler->fb_sequence, 0, sizeof(float)*view_width()*view_height()*3);
  memset(rt.pointsampler->fb_sequence_nc, 0, sizeof(float)*view_width()*view_height()*3);
  for (int i=0;i<NUM_BEAUTY_SPLITS;++i)
  {
    memset(rt.pointsampler->fb_uniform_split[i], 0, sizeof(float)*view_width()*view_height()*3);
    memset(rt.pointsampler->fb_guiding_split[i], 0, sizeof(float)*view_width()*view_height()*3);
  }
  memset(rt.pointsampler->throughput_histogram, 0, sizeof(uint64_t)*rt.pointsampler->num_throughput_bins);
  memset(rt.pointsampler->throughput_histogram_old, 0, sizeof(uint64_t)*rt.pointsampler->num_throughput_bins);
  memset(rt.pointsampler->throughput_histogram_diff, 0, sizeof(uint64_t)*rt.pointsampler->num_throughput_bins);
}

void pointsampler_cleanup(pointsampler_t *s)
{
  guided_cleanup(s->cache);
  free(s->guiding_sampling_hist);
  free(s->guiding_hist);
  free(s->fb_uniform);
  free(s->fb_guiding);
  free(s->fb_sequence);
  free(s->fb_sequence_nc);
  for (int i=0;i<NUM_BEAUTY_SPLITS;++i)
  {
    free(s->fb_uniform_split[i]);
    free(s->fb_guiding_split[i]);
  }
  free(s->throughput_histogram);
  free(s->throughput_histogram_old);
  free(s->throughput_histogram_diff);
  free(s->rand);
  free(s);
}
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}
void pointsampler_finalize(pointsampler_t *s) {}

void pointsampler_stop_learning(pointsampler_t *s)
{
  s->render_mode = rm_render;
#if CLEAR_AFTER_LEARNING==1
  printf("clear frame buffer after learning\n");
  view_clear_frame(rt.view);
  memset(s->fb_sequence, 0, sizeof(float)*view_width()*view_height()*3);
#endif
  guided_build_cdf(s->cache, s->num_nb, s->uniform_ratio, rm_learn);
  guided_build_cdf(s->cache, s->num_nb, s->uniform_ratio, rm_render);
  guided_export_cache_info(s->cache);
}

float pointsampler(path_t *p, int i)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[i];

  int v = p->length;
  const int end = p->v[v].rand_beg;
  const int dim = end + i;
  if(p->index == -1ul || dim >= halton_get_num_dimensions())
    // degenerate to pure random mersenne twister
    return points_rand(rt.points, tid);
  else
    // note that this clips the bits in p->index to 32:
    return halton_sample(&rt.pointsampler->h, dim, p->index);
}

void pointsampler_splat_framebuffer(
    path_t *p, 
    float value,
    sampling_type_t st)
{
#if HEROWAVELENGTH==1
  // actually splat with a couple hero wavelegths
  const int num_wavelengths = 4;
  double pdfsum = 0.0f;
  double val[num_wavelengths];
  double pdf[num_wavelengths];
  double mis[num_wavelengths];
  double lam[num_wavelengths];
  path_t perpath;

  for(int i=0;i<num_wavelengths;i++)
  {
    perpath = *p;
    float l_rand = (p->lambda - spectrum_sample_min)/(spectrum_sample_max - spectrum_sample_min);
    l_rand = fmodf(l_rand + i/(float)num_wavelengths, 1.0f);
    perpath.lambda = lam[i] = spectrum_sample_min + (spectrum_sample_max - spectrum_sample_min)*l_rand;

    for(int v=1; v<p->length; v++)
    {
      if(path_edge_init_volume(&perpath, v)) goto fail;
      // catch volume stack propagation edge case:
      // TODO: make sure these paths aren't passed on until here (inconsistent mode coming from measurement_contrib eval vs sample)
      if(!(perpath.v[v].flags & s_environment) && !(perpath.v[v].mode & s_sensor) &&
          primid_invalid(perpath.v[v].hit.prim) && (perpath.e[v].vol.shader < 0)) goto fail;
      shader_prepare(&perpath, v);
    }
    val[i] = path_measurement_contribution_dwp(&perpath, 0, perpath.length-1);

    const float u = rt.pointsampler->uniform_ratio;
    const double updf = sampler_sum_pdf_dwp(&perpath);
    double gpdf = guided_pdf_dwp(rt.pointsampler->cache, &perpath);
    pdf[i] = u * updf + (1.0-u) * gpdf;
    //pdf[i] = sampler_sum_pdf_dwp(&perpath);
    
    if (st == st_uniform)
      mis[i] = (u * updf) / pdf[i];
    else
      mis[i] = ((1.0-u) * gpdf) / pdf[i];
    //mis[i] = sampler_mis_weight(&perpath);
    pdfsum += pdf[i];
    if(0)
    {
fail:
      pdf[i] = 0.0f;
      val[i] = 0.0f;
    }
  }

  float col[3] = {0.0f};
  for(int i=0;i<num_wavelengths;i++)
  {
    if(!(val[i] > 0.0) || !(pdf[i] > 0.0)) continue;
    double v = val[i] / pdfsum  * mis[i];
    if(!(v == v)) continue;
    float c2[3];
    spectrum_p_to_camera(lam[i], v, c2);
    for(int k=0;k<3;k++) col[k] += c2[k];
  }
  view_splat_col(&perpath, col);
#else
  view_splat(p, value);
#endif
}

void pointsampler_splat(path_t *p, float value)
{
#if 0
  {
  const double mc = path_measurement_contribution_dwp(p, 0, p->length-1);
  const double updf = sampler_sum_pdf_dwp(p);
  if(updf > 0.0) view_splat(p, mc/updf);
  return; // XXX
  }
#endif
  const double mc = path_measurement_contribution_dwp(p, 0, p->length-1);
  
  // this is only called from sampler for uniform samples
  const float u = rt.pointsampler->uniform_ratio;
  const double updf = sampler_sum_pdf_dwp(p);
  double gpdf = guided_pdf_dwp(rt.pointsampler->cache, p);
  const double pdf = u * updf + (1.0-u) * gpdf;

  const double contrib = mc/pdf;
  //printf("t %f, mc %f, pdf %f\n", value, mc, sampler_sum_pdf_dwp(p));
  // fprintf(stderr, "umc %8g %8g %8g => %9g\n", mc, updf, gpdf, contrib);
  if(!(contrib < FLT_MAX)) return;
  
  uint64_t contrib_bin = 0;
  int contrib_scale = 1<<17;
  float contrib_scaled = contrib * contrib_scale;
  if (contrib_scaled >= 1.f)//everything else is in bin 0
  {
    contrib_bin = (uint64_t) log2f(contrib_scaled) + 1;
    if (contrib_bin > rt.pointsampler->num_throughput_bins-1) contrib_bin = rt.pointsampler->num_throughput_bins-1;
  }
  rt.pointsampler->throughput_histogram[contrib_bin]++;

  if (rt.pointsampler->render_mode == rm_learn)
  {
#if STOP_LEARNING==1
    pointsampler_splat_framebuffer(p, contrib, st_uniform);
#endif
    if(p->v[0].mode & s_emit)
    {
      path_t q;
      path_reverse(&q, p);
      guided_record_path(rt.pointsampler->cache, &q, (uint32_t)(-1), mc, gpdf, updf, u, st_uniform);
    }
    else
      guided_record_path(rt.pointsampler->cache, p, (uint32_t)(-1), mc, gpdf, updf, u, st_uniform);
    
  }
#if STOP_LEARNING==1
  else
#endif
  {
#if DUMP_BEAUTY==1
    splat_fb(rt.pointsampler->fb_uniform, p, contrib);
    pointsampler_splat_framebuffer(p, contrib, st_uniform);
#if DUMP_BEAUTY_SPLITS==1
    int split = CLAMP(p->length-3, 0, NUM_BEAUTY_SPLITS-1);
    splat_fb(rt.pointsampler->fb_uniform_split[split], p, contrib);
#endif
#else
    pointsampler_splat_framebuffer(p, contrib, st_uniform);
#endif
  }
  splat_fb(rt.pointsampler->fb_sequence, p, contrib);
  splat_fb(rt.pointsampler->fb_sequence_nc, p, contrib);
  
#if 0
    // visualise pdf without sampling it
    view_splat(p, guided_pdf_dwp(rt.pointsampler->cache, p));
#endif
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  pointsampler_t *s = rt.pointsampler;
  // cast certain percentage of uniform and guided samples
  // at the beginning 100% uniform, decreasing with progressions

  const float u = s->uniform_ratio;
  if(tent->index < s->halton_end)
  { // this is the uniform version:
    path_init(tent, tent->index, tent->sensor.camid);
    tent->index = s->halton_org + (tent->index - s->halton_beg);
    sampler_create_path(tent);
  }
  else
  { // use guided construction:
    path_init(tent, tent->index, tent->sensor.camid);
    tent->index = -1ul; // switch off qmc
    
    // sample a guide path
    guided_samples++;
    int err = 0;
    int pathid = guided_sample_guide_path(s->cache);
    if (pathid < 0)
      return;
    err = guided_sample(s->cache, tent, pathid);
    // if(err>1) fprintf(stderr, "er = %d\n", err);
    if(err <= 1 && tent->throughput > 0)
    {
      // view_splat(tent, 1.0f); // DEBUG see sampling probability
      double gpdf  = guided_pdf_dwp(s->cache, tent);
      //if (gpdf > 0)
      //  gpdf = MAX(1.0, gpdf);
      const double updf = sampler_sum_pdf_dwp(tent);
      const double pdf = u * updf + (1.0-u) * gpdf;
      if(!(gpdf > 0.0)) return; // XXX should be killed during sampling
      // assert(gpdf > 0.0); // we just constructed it, right?
      const double mc = path_measurement_contribution_dwp(tent, 0, tent->length-1);
      const double contrib = mc/pdf;
      // fprintf(stderr, "mc %8g %8g %8g => %9g\n", mc, updf, gpdf, contrib);
      if(!(contrib < FLT_MAX)) return;
      // const double contrib = 1.0/pdf;
      s->guiding_sampling_hist[(int)(tent->sensor.pixel_i)+(int)(tent->sensor.pixel_j)*view_width()]+=1.f;

      uint64_t contrib_bin = 0;
      int contrib_scale = 1<<17;
      float contrib_scaled = contrib * contrib_scale;
      if (contrib_scaled >= 1.f)//everything else is in bin 0
      {
        contrib_bin = (uint64_t) log2f(contrib_scaled) + 1;
        if (contrib_bin > rt.pointsampler->num_throughput_bins-1) contrib_bin = rt.pointsampler->num_throughput_bins-1;
      }
      rt.pointsampler->throughput_histogram[contrib_bin]++;

      if (rt.pointsampler->render_mode == rm_learn)
      {
#if STOP_LEARNING==1
        pointsampler_splat_framebuffer(tent, contrib, st_guided);
#endif
        guided_record_path(s->cache, tent, pathid, mc, gpdf, updf, u, st_guided);
      }
#if STOP_LEARNING==1
      else
#endif
      {
#if DUMP_BEAUTY==1
        splat_fb(rt.pointsampler->fb_guiding, tent, contrib);
        pointsampler_splat_framebuffer(tent, contrib, st_guided);
#else
        pointsampler_splat_framebuffer(tent, contrib, st_guided);
#endif
#if DUMP_BEAUTY_SPLITS==1
        int split = CLAMP(tent->length-3, 0, NUM_BEAUTY_SPLITS-1);
        splat_fb(rt.pointsampler->fb_guiding_split[split], tent, contrib);
#endif
      }
      splat_fb(rt.pointsampler->fb_sequence, tent, contrib);
      splat_fb(rt.pointsampler->fb_sequence_nc, tent, contrib);
    }
    else guided_misses++;
  }
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
	uint64_t histogram_max = 0;
  uint64_t histogram_max_diff = 0;
  uint64_t sum_diff = 0;
  uint64_t sum_hist = 0;
  for(uint64_t bin=0;bin<s->num_throughput_bins;bin++)
  {
    //comment in following line for loglog histogram:
    //s->throughput_histogram[bin] = log2f(s->throughput_histogram[bin]);
    histogram_max = MAX(histogram_max, s->throughput_histogram[bin]);
    sum_hist += s->throughput_histogram[bin];
    if (rt.frames >= 1)
    {
      //Replace old histogram with histogram diff
      uint64_t diff = MAX((double)(s->throughput_histogram_old[bin])/rt.frames, (double)s->throughput_histogram[bin])
                    - MIN((double)(s->throughput_histogram_old[bin])/rt.frames, (double)s->throughput_histogram[bin]);
      s->throughput_histogram_diff[bin] = diff;
      sum_diff += diff;
      histogram_max_diff = MAX(histogram_max_diff, diff);
    }
  }

#if PRINT_THROUPUT_HISTOGRAM == 1
  if (histogram_max <= 0)
  {
    printf("\nhistogram max = %ld\n", histogram_max);
  }
  else
  {
    printf("           ");
    printf("histogram                       ");
    printf("           ");
    printf("histogram avg                   ");
    printf("           ");
    printf("histogram diff (max %09ld)  \n", histogram_max_diff);
    for(int h=lines-1;h>=0;h--)
    {
      printf("           ");
      for(int bin=0;bin<s->num_throughput_bins;bin++)
      {
        float level = s->throughput_histogram[bin] / (float)histogram_max;
        float fill = level*lines-h;
        if(fill <= 0) printf(" ");
        else if(fill <= 1./8.) printf("\u2581");
        else if(fill <= 2./8.) printf("\u2582");
        else if(fill <= 3./8.) printf("\u2583");
        else if(fill <= 4./8.) printf("\u2584");
        else if(fill <= 5./8.) printf("\u2585");
        else if(fill <= 6./8.) printf("\u2586");
        else if(fill <= 7./8.) printf("\u2587");
        else /*if(fill <= 8./8.)*/  printf("\u2588");
      }
      printf("           ");
      for(int bin=0;bin<s->num_throughput_bins;bin++)
      {
        float level = (double)(s->throughput_histogram_old[bin]) / 100000.0 / MAX(rt.frames,1);
        float fill = level*lines-h;
        if (fill <= 0) printf(" ");
        else if(fill <= 1./8.) printf("\u2581");
        else if(fill <= 2./8.) printf("\u2582");
        else if(fill <= 3./8.) printf("\u2583");
        else if(fill <= 4./8.) printf("\u2584");
        else if(fill <= 5./8.) printf("\u2585");
        else if(fill <= 6./8.) printf("\u2586");
        else if(fill <= 7./8.) printf("\u2587");
        else printf("\u2588");
      }
      printf("           ");
      for(int bin=0;bin<s->num_throughput_bins;bin++)
      {
        float level = (double)(s->throughput_histogram_diff[bin]) / 100000.0;
        float fill = level*lines-h;
        if (fill <= 0) printf(" ");
        else if(fill <= 1./8.) printf("\u2581");
        else if(fill <= 2./8.) printf("\u2582");
        else if(fill <= 3./8.) printf("\u2583");
        else if(fill <= 4./8.) printf("\u2584");
        else if(fill <= 5./8.) printf("\u2585");
        else if(fill <= 6./8.) printf("\u2586");
        else if(fill <= 7./8.) printf("\u2587");
        else printf("\u2588");
      }
      printf("\n");
    }
    printf("          1");
    printf("|  5|    |    |  20|    |    |  ");
    printf("          1");
    printf("|  5|    |    |  20|    |    |  ");
    printf("          1");
    printf("|  5|    |    |  20|    |    |  \n");
  }
#endif

  if (rt.frames == 0)
  {
    s->uniform_ratio = 1.f;
    s->learn_time = common_time_wallclock();
    s->start_time = common_time_wallclock();
  }
  else {
    s->uniform_ratio = 0.5f;
  }
#if STOP_LEARNING==1
  double elapsed_time = common_time_wallclock() - s->learn_time;
  //if (guided_num_guide_paths(s->cache) > 2000 && s->render_mode == rm_learn) 
  //if (rt.frames >= LEARN_ITERATIONS && s->render_mode == rm_learn) 
  if (elapsed_time > LEARN_TIME && s->render_mode == rm_learn) 
  {
    pointsampler_stop_learning(s);
    s->learn_time = common_time_wallclock() - s->learn_time;
    s->learn_iterations = rt.frames;
    s->learn_paths = guided_num_guide_paths(s->cache);
    printf("learning finished (%gs, %ld iterations, %ld guide paths)\n", s->learn_time, s->learn_iterations, s->learn_paths);
  }
#endif

  s->num_nb = 10;
  s->heap_size = 100;
  guided_prepare_frame(s->cache, s->num_nb, s->heap_size);

  double time_build_cdf = common_time_wallclock();

  if (s->render_mode == rm_learn)
  {
    guided_build_cdf(s->cache, s->num_nb, s->uniform_ratio, rm_learn);
  }

  time_build_cdf = common_time_wallclock() - time_build_cdf;
  double time_prog = common_time_wallclock() - s->time_stamp;
  s->time_stamp = common_time_wallclock();
  
  if (rt.frames > 1)
    fprintf(stderr, "render_mode %d, num guide paths %d, uniform ratio %f, guided miss : %.02f%%, time prog %g, time build cdf %g (%.2f%%)\n", 
        s->render_mode, guided_num_guide_paths(s->cache), s->uniform_ratio, 
        guided_misses*100.0/guided_samples, time_prog, time_build_cdf, time_build_cdf*100.0/time_prog);

  const uint64_t num = s->uniform_ratio * (rt.threads->end - rt.threads->counter);
  s->halton_org += (s->halton_end - s->halton_beg);
  s->halton_beg = rt.threads->counter;
  s->halton_end = rt.threads->counter + num;

  // would wrap around int limit and run out of bits for radical inverse?  stop
  // to re-init the random bit permutations (we only pass a seed, will be
  // randomised internally):
  if(s->halton_end >> 32 > s->reinit)
    halton_init_random(&s->h, rt.anim_frame + ++s->reinit);
  
  guided_samples = 0;
  guided_misses = 0;

  const uint64_t n = rt.frames - s->learn_iterations;
  if(!(n == 0 || (n & (n-1))))
  {
    //uint32_t time = (uint32_t)(common_time_wallclock() - s->start_time);
    printf("write sequence rt.frames = %lu, n = %lu, learn_it = %lu\n", rt.frames, n, s->learn_iterations);
    {
      char filename[1024];
      snprintf(filename, sizeof(filename), "%s_fb_seq_%06lu", s->render_mode == rm_learn ? "learn" : "render", rt.frames);
      write_fb(s->fb_sequence, filename, view_overlays());
    }
    {
      char filename[1024];
      snprintf(filename, sizeof(filename), "fb_seq_nc_%06lu", rt.frames);
      write_fb(s->fb_sequence_nc, filename, rt.frames);
    }
  }

#if DUMP_BEAUTY==1
  if (s->render_mode == rm_render && (rt.frames % 32 == 0))
  {
    write_fb(s->fb_guiding, "fb_guiding", view_overlays());
    write_fb(s->fb_uniform, "fb_uniform", view_overlays());
#if DUMP_BEAUTY_SPLITS==1
    for (int i=0;i<NUM_BEAUTY_SPLITS;++i)
    {
      {
        char filename[1024];
        snprintf(filename, sizeof(filename), "fb_guiding_split_%03d", i);
        write_fb(s->fb_guiding_split[i], filename, view_overlays());
      }
      {
        char filename[1024];
        snprintf(filename, sizeof(filename), "fb_uniform_split_%03d", i);
        write_fb(s->fb_uniform_split[i], filename, view_overlays());
      }
    }
#endif
  }
#endif
#if DUMP_GUIDING_HIST==1
  if (s->render_mode == rm_learn)
  {
    {
      char filename[1024];
      snprintf(filename, sizeof(filename), "guiding_hist%03ld", rt.frames);
      memset(s->guiding_hist, 0, sizeof(float)*view_width()*view_height()*3);
      guided_gpath_hist(s->cache, s->guiding_hist, view_width(), view_height());
      screenshot_write(filename, s->guiding_hist, 3, 0, 3, view_width(), view_height(), 1.f);
    }
  }
#endif
  for(int bin=0;bin<s->num_throughput_bins;bin++)
    s->throughput_histogram_old[bin] += s->throughput_histogram[bin];
  memset(s->throughput_histogram, 0, sizeof(uint64_t)*s->num_throughput_bins);
}
