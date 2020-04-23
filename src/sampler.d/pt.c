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
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "sampler.h"
#include "spectrum.h"
#include "pointsampler.h"

// std backward pathtracer

typedef struct sampler_t {} sampler_t;
sampler_t *sampler_init() {return 0;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

static inline mf_t sampler_mis(const path_t *p)
{
  // this is just the hero wavelength weight:
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf), mf_set1(mf_hsum(md_2f(pdf))));
}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    if(path->v[path->length-1].mode & s_emit)
    {
      const mf_t w = sampler_mis(path);
      pointsampler_splat(path, mf_mul(w, path_throughput(path)));
      if(path->length > 3)
        if(path_russian_roulette(path, MIN(1.0f, mf(path->v[path->length-1].throughput, 0)/mf(path->v[path->length-2].throughput, 0))))
          return;
    }
  }
}

mf_t sampler_throughput(path_t *path)
{
  if(path->length < 1) return mf_set1(0.0f);
  const md_t measurement = path_measurement_contribution_dx(path, 0, path->length-1);
  if(mf_all(mf_lte(md_2f(measurement), mf_set1(0.0f))))
    return mf_set1(0.0f);
  md_t pdf = md_set1(1.0);
  for(int k=0;k<path->length;k++)
    pdf = md_mul(pdf, mf_2d(path_pdf_extend(path, k)));
  return md_2f(md_div(measurement, pdf));
}

md_t sampler_mis_weight(path_t *p)
{
  return md_set1(1.0);
}

md_t sampler_sum_pdf_dwp(path_t *p)
{
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(mf_div(path_pdf_extend(p, v), mf_set1(path_G(p, v)))));
  return pdf;
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : pathtracer\n");
}
