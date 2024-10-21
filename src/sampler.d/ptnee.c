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
#include "pointsampler.h"
#include "display.h"
#include "shader.h"
#include "pathspace/nee.h"

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
  display_control_add(rt.display, "[ptnee] path verts", &s->max_path_len, 2, PATHSPACE_MAX_VERTS, 1, 0, 1);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    if(path->length == 2 && path->v[2].mode & s_emit)
    { // cannot be created by nee
      mf_t throughput = path_throughput(path);
      if(mf_any(mf_gt(throughput, mf_set1(0.0f))))
        pointsampler_splat(path, throughput);
    }

    if(path->length >= rt.sampler->max_path_len) return;

    if(nee_sample(path)) return;
    const int v2 = path->length-1;
    mf_t throughput = path_throughput(path);
    if(mf_any(mf_gt(throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_emit))
      pointsampler_splat(path, throughput);
    path_pop(path);
  }
}

md_t sampler_mis_weight(path_t *p)
{
  return md_set1(1.0);
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : pathtracer with only next event estimation\n");
}
