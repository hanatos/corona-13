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

#ifndef SAMPLER_LT_H
#define SAMPLER_LT_H

#include "pathspace.h"
#include "pathspace/nee.h"
#include "pathspace/mvnee.h"
#include "pointsampler.h"

// std forward pathtracer (light tracer)

typedef struct sampler_t {} sampler_t;
sampler_t *sampler_init() {return NULL;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

static inline mf_t
sampler_mis(const path_t *p, const mf_t pdf)
{
  // evaluate combined balance heuristic for wavelength and path construction:
  md_t pdf_path = md_set1(1.0);
  for(int v=1;v<p->length-1;v++) // forget about G terms, they'll cancel out:
    pdf_path = md_mul(pdf_path, mf_2d(p->v[v].pdf));

  md_t our = md_mul(mf_2d(pdf), pdf_path);
  return md_2f(md_div(our, md_set1(md_hsum(our))));
}

void sampler_create_path(path_t *path)
{
  path->v[0].mode = s_emit; // start at the light
  while(1)
  {
    if(path_extend(path)) return;
    if(nee_sample(path)) return;
    const int v2 = path->length-1;
    if(mf_all(mf_gt(path->v[v2].throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_sensor))
    {
      const mf_t weight = sampler_mis(path, path->v[v2].pdf);
      pointsampler_splat(path, mf_mul(weight, path->v[v2].throughput));
    }
    path_pop(path);
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : lighttracer\n");
}

#endif
