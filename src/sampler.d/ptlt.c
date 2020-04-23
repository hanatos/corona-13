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
#include "pointsampler.h"

// mis combination of light tracer and path tracer (both with next event estimation)

typedef struct sampler_t {} sampler_t;
sampler_t *sampler_init() {return 0;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

// return mis weight knowing that we only have three techniques here: pt ptdl lt
// light_v is the number of vertices created from the light (i.e. 0, 1, or length-2 respectively)
static inline float sampler_mis(path_t *path, int light_v)
{
  const int light_dir = (path->v[0].mode & s_emit);
  if(path->length == 2)
  {
    if(!light_dir && (light_v == 0)) return 1.0f;
    else return 0.0f; // no other technique is good at sampling directly visible light sources
  }
  // got a path with l=light_v light vertices and e=path->length-light_v eye vertices.

  // get our pdf as sampled:
  double pdf = path_pdf(path);
  // light tracing case, e=1, only able to construct paths len >= 3
  if(light_dir && (light_v == path->length-1))
  {
    // e=0 doesn't exist (can't hit the aperture by chance).
    // start at length-2, because length-1 is already sampled from the camera (lt is the e=1 case).
    // this is special cased in pdf_extend_adjoint.
    double pt_pdf = 1.0;
    for(int k=path->length-2;k>0;k--)
      pt_pdf *= path_pdf_extend_adjoint(path, k);

    const double connect_pdf = path_pdf_next_event_adjoint(path, 0);
    const double extend_pdf  = path_pdf_extend_adjoint(path, 0);

    // weight (power heuristic alpha=2)
    const double mis = pdf*pdf/(pt_pdf*pt_pdf*(connect_pdf*connect_pdf + extend_pdf*extend_pdf) + pdf*pdf);

    return mis;
  }

  double lt_pdf = 1.0f;
  // again, first two vertices from camera go together in one call:
  for(int k=path->length-2;k>0;k--)
    lt_pdf *= path_pdf_extend_adjoint(path, k);

  // last vertex is connection to camera:
  lt_pdf *= path_pdf_next_event_adjoint(path, 0);

  // path tracing case, e=length
  if(light_v == 0)
  {
    // ptdl is almost like ours, only the last vertex didn't scatter but was next event estimation:
    double ptdl_pdf = pdf;
    ptdl_pdf /= path->v[path->length-1].pdf;
    ptdl_pdf *= path_pdf_next_event(path, path->length-1);

    const double mis = pdf*pdf/(lt_pdf*lt_pdf + ptdl_pdf*ptdl_pdf + pdf*pdf);
    return mis;
  }

  // ptdl case, e=length-1
  if(light_v == 1)
  {
    // pt pdf used extend on the last vertex, not next event estimation
    double pt_pdf = pdf;
    pt_pdf /= path->v[path->length-1].pdf;
    pt_pdf *= path_pdf_extend(path, path->length-1);

    const double mis = pdf*pdf/(lt_pdf*lt_pdf + pt_pdf*pt_pdf + pdf*pdf);
    return mis;
  }

  // that was a weird technique. shouldn't be in this file.
  return 0.0f;
}

void sampler_create_path(path_t *path)
{
  // lt path
  path->v[0].mode = s_emit; // start at the light
  while(1)
  {
    if(path_extend(path)) break;
    if(path_next_event(path)) break;
    const int v2 = path->length-1;
    if(path->v[v2].throughput > 0.0f && (path->v[v2].mode & s_sensor))
    {
      const float weight = sampler_mis(path, path->length-1);
      pointsampler_splat(path, path->v[v2].throughput * weight);
    }
    path_pop(path);
  }

  // eye path
  path_init(path);
  // FIXME: restore v[0].rand_beg to work with kelemen mlt? or randomly decide lt vs pt!
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
      const float weight = sampler_mis(path, 0);
      pointsampler_splat(path, path_throughput(path) * weight);
    }
    if(path->v[v].flags & s_environment) return;
    if(path->v[v].flags & s_emit) return;

    if(path_next_event(path)) break;
    const int v2 = path->length-1;
    if(path->v[v2].throughput > 0.0f && (path->v[v2].mode & s_emit))
    {
      const float weight = sampler_mis(path, 1);
      pointsampler_splat(path, path->v[v2].throughput * weight);
    }
    path_pop(path);
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : mis-combined pathtracer with next event estimation and lighttracer\n");
}
