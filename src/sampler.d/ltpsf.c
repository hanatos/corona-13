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
#include "pathspace/equiangular.h"
#include "pointsampler.h"

// std forward pathtracer (light tracer)

// #define DEBUG_COMP_MVNEE=1
// #define DEBUG_COMP_EQUI=1

typedef struct sampler_t {} sampler_t;
sampler_t *sampler_init() {return NULL;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

// compute HWL MIS weight
static inline mf_t
sampler_mis(
    const path_t *p,
    const int     tech) // 0 - this path fwd scattering, 1 - nee, 2 - mvnee
{
  // compute product of vertex area pdfs along the path excluding the last two vertices ve and ve-1:
  const int ve = p->length - 1;
  md_t pdf_path = md_set1(1.0);
  for(int v=1;v<ve-1;v++) // the first two vertices 0,1 are created in a single path_extend call and have only 1 pdf
    pdf_path = md_mul(pdf_path, mf_2d(p->v[v].pdf));

  md_t pdf1 = md_mul(mf_2d(path_pdf_extend(p, ve-1)), mf_2d(nee_pdf(p, ve)));
#ifdef DEBUG_COMP_MVNEE
  md_t pdf2 = mf_2d(mvnee_pdf(p, ve));
#else
  md_t pdf2 = md_set1(0.0);
#endif
#ifdef DEBUG_COMP_EQUI
  md_t pdf3 = mf_2d(equiangular_pdf(p, ve));
#else
  md_t pdf3 = md_set1(0.0);
#endif
  md_t our;

#if 0 // trust the pdf on the path
  our = md_mul(mf_2d(p->v[ve-1].pdf), mf_2d(p->v[ve].pdf));
  if(tech == 1) pdf1 = our;// = pdf1;
  if(tech == 2) pdf2 = our;// = pdf2;
  if(tech == 3) pdf2 = our;// = pdf3;
#else // trust the explicit eval
  if(tech == 1) our = pdf1;
  if(tech == 2) our = pdf2;
  if(tech == 3) our = pdf3;
#endif

#if 0 // "power heuristic"
  pdf1 = md_mul(pdf1, pdf1);
  pdf2 = md_mul(pdf2, pdf2);
  pdf3 = md_mul(pdf3, pdf3);
  our  = md_mul(our, our);
#endif
  // evaluate combined balance heuristic for wavelength and path construction:
  return mf_div(md_2f(our), mf_set1(mf_hsum(md_2f(md_add(pdf1, md_add(pdf2, pdf3))))));

#if 0
  // our path has the following two vertex area pdfs at [ve-1] and [ve]:
  md_t our = md_mul(md_mul(mf_2d(p->v[ve-1].pdf), mf_2d(p->v[ve].pdf)), pdf_path);
  md_t other;
  if(tech == 1)                         // we are nee, the other is mvnee
    other = mf_2d(mvnee_pdf(p, ve)); // this returns the product of the last two vertices
  else if(tech == 2)                    // we are mvnee, the other is nee
    other = md_mul(mf_2d(p->v[ve-1].pdf), mf_2d(nee_pdf(p, ve)));
  else return mf_set1(0.0f);
  // multiply to rest of path:
  other = md_mul(pdf_path, other);

  // evaluate combined balance heuristic for wavelength and path construction:
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
// #define WHALE
// #define REF
// #ifdef REF
      if(nee_sample(path)) return;
      const int v2 = path->length-1;
// #ifdef WHALE // whale:
//       if(path->length == 5)                   // 5-vtx paths with:
//       // if(path->v[2].hit.prim.shapeid == 0)    // caustic on whale
//       if(primid_invalid(path->v[3].hit.prim)) // scattered once in the medium
//       // if(dotproduct(path->e[3].omega, path->e[4].omega) >= 0.0) // only account for forward scattering paths (as mvnee)
// #else // sphere:
//       if(path->length == 4)                   // 4-vtx paths with:
//       // if(path->v[1].hit.prim.shapeid == 0)    // highlight on one of the spheres
//       if(primid_invalid(path->v[2].hit.prim)) // scattered once in the medium
//       // if(dotproduct(path->e[2].omega, path->e[3].omega) >= 0.0) // only account for forward scattering paths (as mvnee)
// #endif
      if(mf_any(mf_gt(path->v[v2].throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_sensor))
      {
        mf_t w = sampler_mis(path, 1);
        pointsampler_splat(path, mf_mul(w, path->v[v2].throughput));
      }
      path_pop(path);
// #endif
    }
// #ifndef REF // mvnee
// #ifdef WHALE // whale
//       if(path->length == 3)                   // 3-vtx path prefix (whale)
//       // if(path->v[2].hit.prim.shapeid == 0)    // forming caustic on whale
// #else // spheres
//       if(path->length == 2)                   // 2-vtx path prefix (sphere)
//       // if(path->v[1].hit.prim.shapeid == 0)    // forming highlight on sphere
// #endif
#ifdef DEBUG_COMP_MVNEE
    if(mvnee_possible(path, path->length-1))
    { // mvnee
      if(mvnee_sample(path)) return;
      const int v3 = path->length-1;
      // if(path->v[2].hit.prim.shapeid == 0)
      if(mf_any(mf_gt(path->v[v3].throughput, mf_set1(0.0f))) && (path->v[v3].mode & s_sensor))
      {
        mf_t w = sampler_mis(path, 2);
        pointsampler_splat(path, mf_mul(path->v[v3].throughput, w));
      }
      path_pop(path);
    }
#endif
#ifdef DEBUG_COMP_EQUI
    if(equiangular_possible(path, path->length-1))
    { // equiangular sampling
      if(equiangular_sample(path)) return;
      const int v3 = path->length-1;
      // if(path->v[2].hit.prim.shapeid == 0)
      if(mf_any(mf_gt(path->v[v3].throughput, mf_set1(0.0f))) && (path->v[v3].mode & s_sensor))
      {
        mf_t w = sampler_mis(path, 3);
        pointsampler_splat(path, mf_mul(path->v[v3].throughput, w));
      }
      path_pop(path);
    }
#endif
// #ifdef WHALE // whale
//     if(path->length == 3) return; // done everything we could
// #else // spheres
//     if(path->length == 2) return; // done everything we could
// #endif
// #endif
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : psf lighttracer\n");
}

#endif
