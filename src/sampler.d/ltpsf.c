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

// compute HWL MIS weight
static inline mf_t sampler_mis(const path_t *p, const mf_t pdf, const mf_t pdf2)
{
  md_t pdf_path = md_set1(1.0);
  for(int v=1;v<p->length;v++) // forget about G terms, they'll cancel out:
    pdf_path = md_mul(pdf_path, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf_path), mf_set1(mf_hsum(md_2f(pdf_path))));

  // XXX TODO: for long enough paths, compute until length-2 and eval extend/nee + mvnee for the last two vertices!
  // XXX TODO: need to call from outside with v.pdf * v.pdf and mvnee_pdf or path_pdf_extend * nee_pdf
#if 0 // TODO: fix this up for NEE + mvnee in case we need it
  // evaluate combined balance heuristic for wavelength and path construction:
  md_t pdf_path = md_set1(1.0);
  for(int v=1;v<p->length-1;v++) // forget about G terms, they'll cancel out:
    pdf_path = md_mul(pdf_path, mf_2d(p->v[v].pdf));

  md_t our   = md_mul(mf_2d(pdf), pdf_path);
  md_t other = md_mul(mf_2d(pdf2), pdf_path);
  return mf_div(md_2f(our), mf_set1(mf_hsum(md_2f(md_add(other, our)))));
#endif
}

void sampler_create_path(path_t *path)
{
  path->v[0].mode = s_emit; // start at the light
  while(1)
  {
    if(path_extend(path)) return;

    // only do next event estimation if the anticipated path couldn't have been constructed via mvnee:
    // if(!mvnee_possible(path, path->length-2) || !primid_invalid(path->v[path->length-1].hit.prim))
    { // nee reference
#define WHALE
// #define REF
#ifdef REF
      if(nee_sample(path)) return;
      const int v2 = path->length-1;
#ifdef WHALE // whale:
      if(path->length == 5)                   // 5-vtx paths with:
      if(path->v[2].hit.prim.shapeid == 0)    // caustic on whale
      if(primid_invalid(path->v[3].hit.prim)) // scattered once in the medium
      if(dotproduct(path->e[3].omega, path->e[4].omega) >= 0.0) // only account for forward scattering paths (as mvnee)
#else // sphere:
      if(path->length == 4)                   // 4-vtx paths with:
      if(path->v[1].hit.prim.shapeid == 0)    // highlight on one of the spheres
      if(primid_invalid(path->v[2].hit.prim)) // scattered once in the medium
      if(dotproduct(path->e[2].omega, path->e[3].omega) >= 0.0) // only account for forward scattering paths (as mvnee)
#endif
      if(mf_all(mf_gt(path->v[v2].throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_sensor))
      {
        mf_t w = sampler_mis(path, mf_set1(0.0f), mf_set1(0.0f));
        pointsampler_splat(path, mf_mul(w, path->v[v2].throughput));
      }
      path_pop(path);
#endif
    }
#ifndef REF // mvnee
#ifdef WHALE // whale
      if(path->length == 3)                   // 3-vtx path prefix (whale)
      if(path->v[2].hit.prim.shapeid == 0)    // forming caustic on whale
#else // spheres
      if(path->length == 2)                   // 2-vtx path prefix (sphere)
      if(path->v[1].hit.prim.shapeid == 0)    // forming highlight on sphere
#endif
    // if(mvnee_possible(path, path->length-1))
    { // mvnee
      if(mvnee_sample(path)) return;
      const int v3 = path->length-1;
      // if(path->v[2].hit.prim.shapeid == 0)
      if(mf_all(mf_gt(path->v[v3].throughput, mf_set1(0.0f))) && (path->v[v3].mode & s_sensor))
      {
        mf_t w = sampler_mis(path, mf_set1(0.0f), mf_set1(0.0f));
        pointsampler_splat(path, mf_mul(path->v[v3].throughput, w));
      }
      path_pop(path);
    }
#ifdef WHALE // whale
    if(path->length == 3) return; // done everything we could
#else // spheres
    if(path->length == 2) return; // done everything we could
#endif
#endif
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : psf lighttracer\n");
}

#endif
