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

#ifndef SAMPLER_PTDL1_H
#define SAMPLER_PTDL1_H

#include "pathspace.h"
#include "pointsampler.h"
#include "shader.h"
#include <assert.h>

// std backward pathtracer with next event estimation
// only return one path for metropolis, use rejection sampling:
//
// XXX i cannot make this work with mlt:
// call extend every time
// if hit a light:
//   accum with mis
//   rr to terminate path
// use roughness for roussian roulette to sample nee or not
// if nee:
//   this path ends here, either black or
//   with light source contribution (mis, adjusted by rr weight)
// else: adjust throughput by weight
// XXX so for now it just does next event estimation and /no/ extend connection.
// this is kind of a good fit for glossy only methods (hmc). :(

typedef struct sampler_t {} sampler_t;
sampler_t *sampler_init() {return NULL;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

float sampler_mis(const float pdf, const float pdf2)
{
  if(pdf2 < 0.0f && pdf >= 0.0f) return 0.0f; // other technique is a multiple of a dirac
  if(pdf < 0.0f && pdf2 >= 0.0f) return 1.0f; // we are a multiple of a dirac
  return pdf*pdf/(pdf*pdf + pdf2*pdf2);
}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
      // XXX const float weight = sampler_mis(path->v[v].pdf, path_pdf_next_event(path, v));
      // XXX pointsampler_splat(path, path_throughput(path) * weight);
      return;

      // FIXME: this doesn't work with reflecting light sources/emitting media!
      // if(path->length > 3)
        // if(path_russian_roulette(path, fminf(1.0, path->v[path->length-1].throughput/path->v[path->length-2].throughput)))
          // return;
    }

    const float p_nee = fminf(.5f, path->v[path->length-1].shading.roughness);
    // XXX duplicate use of this random number? technically belongs to v[l-1].
    if(pointsampler(path, s_dim_russian_r) < p_nee)
    {
      if(path_next_event(path)) return;
      const int v2 = path->length-1;
      if(path->v[v2].throughput > 0.0f && (path->v[v2].mode & s_emit))
      {
        const float weight = 1.0;//sampler_mis(path->v[v2].pdf, path_pdf_extend(path, v2));
        pointsampler_splat(path, path->v[v2].throughput * weight/p_nee);
      }
      return; // can only return one path, not a group, so stop here.
    }
    path->v[path->length-1].throughput *= 1.0f / (1.0f-p_nee);
  }
}

float sampler_throughput(path_t *path)
{
  if(path->length < 2) return 0.0f;
  assert(path->v[0].mode & s_sensor);
  float throughput = path_throughput_extend(path, 0) * path_throughput_extend(path, 1);
  for(int k=2;k<path->length-1;k++)
    throughput *= path_throughput_extend(path, k) / (1.0 - fminf(.5f, path->v[k].shading.roughness));

  // FIXME: this doesn't work with reflecting light sources/emitting media!

  const int v = path->length-1;
  if(path->length == 2)
    return throughput * path_throughput_extend(path, v) * path_eval_end(path, v);
  // have to create vertex v still:
  // XXX for participating media with hereogeneous data, it may pay off to call throughput once and check for perfect importance sampling.
  // XXX now we might be stepping through the same voxels twice for transmittance and pdf (and get equal numbers up to mu_t)
  const float measurement = shader_brdf(path, v-1) * shader_vol_transmittance(path, v) * path_eval_end(path, v) * path_G(path, v);
  // const float pdf_extend = shader_pdf(path, v-1) * shader_vol_pdf(path, v) * path_G(path, v);
  const float p_nee = fminf(.5f, path->v[v-1].shading.roughness);
  const float pdf_nee = p_nee * path_pdf_next_event(path, v);
  // XXX where is the MIS weight?
  // return throughput * measurement / ((1.0f-p_nee) * pdf_extend + p_nee * pdf_nee);
  return throughput * measurement / pdf_nee;
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : single-connection pathtracer with next event estimation and mis\n");
}

#endif
