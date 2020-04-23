/*
    This file is part of corona-13.

    copyright (c) 2016 johannes hanika

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

#include "sampler.h"
#include "points.h"
#include "pointsampler.h"
#include "display.h"
#include "shader.h"
#include "spectrum.h"
#include "dbor.h"
#include "pathspace/nee.h"
#include "pathspace/tech.h"

#define FILTER_FIREFLYS 1
#define TRUST_THR 0.25

// std backward pathtracer with next event estimation
static void splat_fb(float* fb, path_t* p, float val)
{
  float col[3];
  spectrum_p_to_camera(p->lambda, val, col);
  const uint64_t offset = (uint64_t)p->sensor.pixel_i+(uint64_t)p->sensor.pixel_j*view_width();
  for(int k=0;k<3;k++) common_atomic_add(fb+3*offset+k, col[k]);
}

static void write_fb(float* fb, const char* filename)
{
  const float inv_o = 1.f/view_overlays();
  screenshot_write(filename, fb, 3, 0, 3, view_width(), view_height(), inv_o);
}

typedef struct sampler_t
{
  float max_path_len;
  dbor_t *dbor;
  
  float *fb_all;
  float *fb_result;
  float *fb_filtered;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  s->max_path_len = PATHSPACE_MAX_VERTS;
  s->dbor = dbor_init(view_width(), view_height(), 1, 20);
  s->fb_filtered = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_result = calloc(view_width()*view_height()*3, sizeof(float));
  s->fb_all = calloc(view_width()*view_height()*3, sizeof(float));
  display_control_add(rt.display, "[ptdl] path verts", &s->max_path_len, 2, PATHSPACE_MAX_VERTS, 1, 0, 1);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  dbor_cleanup(s->dbor);
  free(s->fb_all);
  free(s->fb_result);
  free(s->fb_filtered);
  free(s);
}

void sampler_prepare_frame(sampler_t *s) 
{
  /*
#if FILTER_FIREFLYS==1
  if (rt.frames == 5)
  {
    dbor_export(s->dbor, "ptdl_dbor", view_overlays()+1);
    sampler_clear(s);
    view_clear_frame(rt.view);
  }
#endif
  */
  write_fb(s->fb_filtered, "fb_filtered");
  write_fb(s->fb_all, "fb_all");
  write_fb(s->fb_result, "fb_result");
  if (rt.frames == 10)
    dbor_export(s->dbor, "ptdl_dbor", view_overlays()+1);
}

void sampler_clear(sampler_t *s) {
  memset(s->fb_filtered, 0, sizeof(float)*view_width()*view_height()*3);
  memset(s->fb_result, 0, sizeof(float)*view_width()*view_height()*3);
  memset(s->fb_all, 0, sizeof(float)*view_width()*view_height()*3);
}

double sampler_pdf(path_t *p)
{
  double pdf = 1.0;
  for(int v=0;v<p->length-1;v++)
    pdf *= path_pdf_extend(p, v);
  double pdf_extend = path_pdf_extend(p, p->length-1);
  double pdf_nee    = nee_pdf(p, p->length-1);
  return pdf * (pdf_extend + pdf_nee);
}

double sampler_mis_weight(path_t *p)
{
  double our_pdf = 0.0;
  double pdf_extend = path_pdf_extend(p, p->length-1);
  double pdf_nee    = nee_pdf(p, p->length-1);
  if(p->v[p->length-1].tech == s_tech_nee)
    our_pdf = pdf_nee;
  else if(p->v[p->length-1].tech == s_tech_extend)
    our_pdf = pdf_extend;

  return our_pdf*our_pdf / (pdf_nee*pdf_nee + pdf_extend*pdf_extend);
}

static inline float sampler_mis(const float pdf, const float pdf2)
{
  return pdf*pdf/(pdf*pdf + pdf2*pdf2);
}

static inline float nee_probability(const int v)
{
  return 1;// switch off sub sampling of nee
  if(v <= 2) return 1;
  return 0.1;
}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
      const float weight = sampler_mis(path->v[v].pdf, nee_probability(v-1)*nee_pdf(path, v));
      const float contrib = path_throughput(path) * weight;
      if (contrib > 0.0f)
      {
#if FILTER_FIREFLYS==1
        int firefly = 0;
        if (contrib >= 8.f)
        {
          const float trust = dbor_splat(rt.sampler->dbor, path->sensor.pixel_i, path->sensor.pixel_j, contrib);
          const float thr = MAX(TRUST_THR, 4.f/(rt.frames+1));
          firefly = !(trust > thr * (rt.frames+1));
        }
#else
        int firefly = 0;
#endif
        if (firefly)
        {
          splat_fb(rt.sampler->fb_filtered, path, contrib);
        }
        else
        {
          splat_fb(rt.sampler->fb_result, path, contrib);
          pointsampler_splat(path, contrib);
        }
        splat_fb(rt.sampler->fb_all, path, contrib);
      }
    }
    if(path->length >= rt.sampler->max_path_len) return;

    const float rr = nee_probability(path->length);
    if(points_rand(rt.points, common_get_threadid()) < rr)
    {
    if(nee_sample(path)) return;
    const int v2 = path->length-1;
    float throughput = path_throughput(path) / rr;
    if(throughput > 0.0f && (path->v[v2].mode & s_emit))
    {
      const float weight = sampler_mis(rr * path->v[v2].pdf, path_pdf_extend(path, v2));
      const float contrib = throughput * weight;
      if (contrib > 0)
      {
#if FILTER_FIREFLYS==1
        int firefly = 0;
        if (contrib >= 8.f)
        {
          const float trust = dbor_splat(rt.sampler->dbor, path->sensor.pixel_i, path->sensor.pixel_j, contrib);
          const float thr = MAX(TRUST_THR, 4.f/(rt.frames+1));
          firefly = !(trust > thr * (rt.frames+1));
        }
#else
        int firefly = 0;
#endif
        if (firefly)
        {
          splat_fb(rt.sampler->fb_filtered, path, contrib);
        }
        else
        {
          splat_fb(rt.sampler->fb_result, path, contrib);
          pointsampler_splat(path, contrib);
        }
        splat_fb(rt.sampler->fb_all, path, contrib);
      }
    }
    path_pop(path);
    }
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : pathtracer with next event estimation and mis\n");
#ifdef FNEE
  fprintf(fd, "           using forward next event\n");
#endif
#ifdef SEGMENT_EMISSION
  fprintf(fd, "           using segment emission\n");
#endif
}
