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
#include "points.h"
#include "ext/halton/halton.h"
#include "pathspace/tech.h"
#include "spectrum.h"
#include "shader.h"

#include <stdio.h>
#include <float.h>

typedef struct fake_randoms_t
{
  int enabled;
  float rand[400];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
  halton_t h;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: halton points and hero wavelengths\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
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
  int v = p->length;
  const int end = p->v[v].rand_beg;
  const int dim = end + i;
  if(dim >= halton_get_num_dimensions())
    // degenerate to pure random mersenne twister
    return points_rand(rt.points, common_get_threadid());
  else
    return halton_sample(&rt.pointsampler->h, dim, p->index);
}

void pointsampler_splat(path_t *p, float value)
{
  // actually splat with a couple hero wavelegths
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
  path_t perpath;

  for(int i=0;i<num_wavelengths;i++)
  {
    perpath = *p;
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
      shader_prepare(&perpath, v);
    }
    val[i] = path_measurement_contribution_dwp(&perpath, 0, perpath.length-1);
    pdf[i] = sampler_sum_pdf_dwp(&perpath);
    // val[i] = path_measurement_contribution_dx(&perpath, 0, perpath.length-1);
    // pdf[i] = path_tech_pdf_as_sampled(&perpath);
    mis[i] = sampler_mis_weight(&perpath);
    pdfsum += pdf[i];
    if(0)
    {
fail:
      pdf[i] = 0.0f;
      val[i] = 0.0f;
    }
  }

  for(int i=0;i<num_wavelengths;i++)
  {
    perpath.lambda = lam[i];

    if(!(val[i] > 0.0) || !(pdf[i] > 0.0)) continue;
    double v = val[i] / pdfsum  * mis[i];
    render_splat(&perpath, v);
  }
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  path_init(tent, tent->index, tent->sensor.camid);
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

void pointsampler_prepare_frame(pointsampler_t *s) {}
void pointsampler_stop_learning(pointsampler_t *s) {}
