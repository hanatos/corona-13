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
#include "pathspace/nee.h"
#include "pathspace/tech.h"

// std backward pathtracer with next event estimation

typedef struct sampler_t
{
  float max_path_len;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  s->max_path_len = PATHSPACE_MAX_VERTS;
  display_control_add(rt.display, "[ptdl] path verts", &s->max_path_len, 2, PATHSPACE_MAX_VERTS, 1, 0, 1);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

md_t sampler_sum_pdf_dwp(path_t *p)
{
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length-1;v++)
    pdf = md_mul(pdf, md_div(mf_2d(path_pdf_extend(p, v)), md_set1(path_G(p, v))));
  md_t G = md_set1(path_G(p, p->length-1));
  md_t pdf_extend = mf_2d(path_pdf_extend(p, p->length-1));
  md_t pdf_nee    = mf_2d(nee_pdf(p, p->length-1));
  md_t pdf_other  = md_div(md_add(pdf_extend, pdf_nee), G);
  return md_mul(pdf, pdf_other);
}

md_t sampler_mis_weight(path_t *p)
{
  md_t our_pdf = md_set1(0.0);
  md_t pdf_extend = mf_2d(path_pdf_extend(p, p->length-1));
  md_t pdf_nee    = mf_2d(nee_pdf(p, p->length-1));
  if(p->v[p->length-1].tech == s_tech_nee)
    our_pdf = pdf_nee;
  else if(p->v[p->length-1].tech == s_tech_extend)
    our_pdf = pdf_extend;

  return md_div(md_mul(our_pdf, our_pdf), md_add(md_mul(pdf_nee, pdf_nee), md_mul(pdf_extend, pdf_extend)));
}

static inline mf_t sampler_mis(const path_t *p, const mf_t pdf, const mf_t pdf2)
{
  // evaluate combined balance heuristic for wavelength and path construction:
  md_t pdf_path = md_set1(1.0);
  for(int v=1;v<p->length-1;v++) // forget about G terms, they'll cancel out:
    pdf_path = md_mul(pdf_path, mf_2d(p->v[v].pdf));

  md_t our   = md_mul(mf_2d(pdf), pdf_path);
  md_t other = md_mul(mf_2d(pdf2), pdf_path);
  return mf_div(md_2f(our), mf_set1(mf_hsum(md_2f(md_add(other, our)))));
}

static inline float nee_probability(const path_t *p, const int v)
{
  return 1;// switch off sub sampling of nee

  if(v <= 2) return 1; // protect direct illumination + first bounce

#if 0
  // super hacky way to turn off nee in the inner core
  if(primid_invalid(p->v[v].hit.prim)) //&& (fabsf(p->v[v].interior.mean_cos)<0.1))
  {
    // as soon as the inner volume with mean cos = 0 overlaps at all, don't do nee:
    for(int k=0;k<p->v[v].interior.num_lobes;k++)
      if(p->v[v].interior.l_g[k] == 0.0)
        return 0;
  }
  return 1;
#endif

  // only 10% of nee events:
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
      const mf_t weight = sampler_mis(path, path->v[v].pdf, mf_mul(mf_set1(nee_probability(path, v-1)), nee_pdf(path, v)));
      pointsampler_splat(path, mf_mul(path_throughput(path), weight));
      // if(path->length > 3)
      //   if(path_russian_roulette(path, MIN(1.0, path->v[path->length-1].throughput/path->v[path->length-2].throughput)))
      //     return;
    }
    if(path->length >= rt.sampler->max_path_len) return;

#if 0
    if(path->length > 5)
      if(path_russian_roulette(path, MIN(1.0, path->v[path->length-1].throughput/path->v[path->length-2].throughput)))
        return;
#endif
    // if(path_russian_roulette(path, 0.95f))
      // return;

    const float rr = nee_probability(path, v);
    if(points_rand(rt.points, common_get_threadid()) < rr)
    {
      if(nee_sample(path)) return;
      const int v2 = path->length-1;
      mf_t throughput = mf_div(path_throughput(path), mf_set1(rr));
      if(mf_any(mf_gt(throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_emit))
      {
        const mf_t weight = sampler_mis(path, mf_mul(mf_set1(rr), path->v[v2].pdf), path_pdf_extend(path, v2));
        pointsampler_splat(path, mf_mul(throughput, weight));
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
