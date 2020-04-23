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
#include "ext/halton/halton.h"
#include "pathspace.h"
#include "pathspace/manifold.h"
#include "pathspace/raydifferentials.h"
#include "pathspace/tech.h"
#include "points.h"
#include "pointsampler.h"
#include "render.h"
#include "sampler.h"
#include "shader.h"
#include "spectrum.h"
#include "threads.h"
#include "view.h"

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
  uint64_t reinit;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: halton points, hero wavelengths, and half vector splats, oh my!\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  halton_init_random(&s->h, frame);
  s->reinit = 0;
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
    path_copy(&perpath, p);
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
    mis[i] = sampler_mis_weight(&perpath);
    pdfsum += pdf[i];
    if(0)
    {
fail:
      pdf[i] = 0.0f;
      val[i] = 0.0f;
    }
  }


  if (p->length < 3)
  { // splat as usual and get out
    float col[3] = {0.0f};
    for(int i=0;i<num_wavelengths;i++)
    {
      if(!(val[i] > 0.0) || !(pdf[i] > 0.0)) continue;
      double v = val[i] / pdfsum  * mis[i];
      if(!(v == v)) continue;
      float c2[3];
      spectrum_p_to_camera(lam[i], v, c2);
      for(int k=0;k<3;k++) col[k] += c2[k];
    }
    view_splat_col(&perpath, col);
    return;
  }

  // now determine splatting gaussian for current path (ignore wavelength)
  // we have a path cam -> sp -> emitter now perturb halfvector and compute
  // splatting kernel
  float R[9*PATHSPACE_MAX_VERTS];
  const int s = 0;
  const int e = perpath.length-1;
  perpath.v[s].diffgeo.type = s_pinned_position;
  // TODO: envmaps?
  perpath.v[e].diffgeo.type = s_pinned_position;
  for (int i = s+1; i < e; ++i)
    perpath.v[i].diffgeo.type = s_free;

  if (manifold_compute_tangents(&perpath, s, e) == 0)
    return;

  float mu[2] = {0.0f}, S[4] = {0.0f};

#if 0 // moving endpoint chain
  int err;
  if((err = raydifferentials_vx(&perpath, e, 1, R)))
    return;
  // now R is a 2x2 matrix mapping pixel offsets to vertex dislocation at v[e].
  // multiply that by light source size:
  const float size = 0.8f; // TODO: get from more useful place
  float rd[4] = {
    size * R[0], size * R[1],
    size * R[2], size * R[3],
  };
  mat2_invert(rd, S);

#else // fixed endpoint chain
  double dh_dx = raydifferentials_compute_rd_h(&perpath, R, s, e);
  if(dh_dx == 0.0)
    return;

  float rd_i[3], rd_j[3];
  if(raydifferentials_v1(&perpath, 1.0, 1.0, rd_i, rd_j))
    return;

  // express rd_i and rd_j not in world space but in tangent space of v1 (actually only 2d now):
  const float rduv_i[2] = {dotproduct(rd_i, perpath.v[1].diffgeo.dpdu), dotproduct(rd_i, perpath.v[1].diffgeo.dpdv)};
  const float rduv_j[2] = {dotproduct(rd_j, perpath.v[1].diffgeo.dpdu), dotproduct(rd_j, perpath.v[1].diffgeo.dpdv)};

  const float DsDx[4] = { rduv_i[0], rduv_j[0],
                          rduv_i[1], rduv_j[1] };
  float DxDs[4];
  mat2_invert(DsDx, DxDs);

  // determine splatting kernel
  for (int i = s+1; i < e; ++i)
  {
    float DhDx_i[4];
    float Rd[4] = {R[9*i], R[9*i+1], R[9*i+3], R[9*i+4]};
    mat2_invert(Rd, DhDx_i);

    // what about anisotropic bsdfs?
    float dh[2] = { -perpath.v[i].diffgeo.h[0], -perpath.v[i].diffgeo.h[1] }; // current position of half vector

    float step = perpath.v[i].shading.roughness;
    float dh0[2] = { step, 0.f }; // frame of beckmann bsdf
    float dh1[2] = { 0.f, step };

    // convert current h to tangent frame of v[1] and then screen space ds
    float ds[2], dx[2];
    mat2_mulv(DhDx_i, dh, dx);
    mat2_mulv(DxDs, dx, ds);

    // convert one half vec axis to screen space ds0
    float ds0[2], dx0[2];
    mat2_mulv(DhDx_i, dh0, dx0);
    mat2_mulv(DxDs, dx0, ds0);

    // other half vec axis to screen space ds1
    float ds1[2], dx1[2];
    mat2_mulv(DhDx_i, dh1, dx1);
    mat2_mulv(DxDs, dx1, ds1);

    // create gaussian parameters:
    float mu_i[2] = {ds[0], ds[1]};
    float S_i[4] = {
      ds0[0]*ds0[0] + ds1[0]*ds1[0],
      ds0[0]*ds0[1] + ds1[0]*ds1[1],
      ds0[0]*ds0[1] + ds1[0]*ds1[1],
      ds0[1]*ds0[1] + ds1[1]*ds1[1]
    };

    // convolve with what we already have:
    for(int k=0;k<2;k++) mu[k] += mu_i[k];
    for(int k=0;k<4;k++) S [k] += S_i [k];
  }

  // TODO: what about light source, is it a gaussian, too?
  // TODO: what about brightness, do we really just want to normalise it?
#endif

#if 1
  // multiply by large envelope gaussian around pixel:
  const float sigma = 2.0f;
  float S_p[4] = {sigma*sigma, 0.0f, 0.0f, sigma*sigma};
  float Si_p[4], Si[4];
  mat2_invert(S, Si);
  mat2_invert(S_p, Si_p);
  float St[4], mu_t[2];
  for(int k=0;k<4;k++) St[k] = Si[k] + Si_p[k];
  // overwrite S:
  mat2_invert(St, S);
  // now the mean:
  mat2_mulv(Si, mu, mu_t);
  mat2_mulv(S, mu_t, mu);
  // float mu_p[2] = {0.0f};
  // mat2_mulv(Si_p, mu_p, mu_t); // this is currently all 0
  // mat2_mulv(S, mu_t, mu_p);
  // for(int k=0;k<2;k++) mu[k] += mu_p[k]; // and this not needed
#endif
  
#if 1
  const float t1 = view_sample_time(p, 1.0f);
  // retime first intersection to derive motion
  assert(p->v[0].mode & s_sensor); // didn't implement for light tracing yet
  float i0, j0, i1, j1;
  path_t rdpath;
  path_copy(&rdpath, p);
  rdpath.sensor = p->sensor;
  // shutter open
  prims_retime(rt.prims, p->v[1].hit.prim, &rdpath.v[1].hit, 0.0f);
  rdpath.time = 0.0f;
  view_cam_connect(&rdpath);
  i0 = rdpath.sensor.pixel_i - p->sensor.pixel_i;
  j0 = rdpath.sensor.pixel_j - p->sensor.pixel_j;
  // shutter close
  prims_retime(rt.prims, p->v[1].hit.prim, &rdpath.v[1].hit, t1);
  rdpath.time = 1.0f;
  view_cam_connect(&rdpath);
  i1 = rdpath.sensor.pixel_i - p->sensor.pixel_i;
  j1 = rdpath.sensor.pixel_j - p->sensor.pixel_j;
  // fit a gaussian

  float mu_m[2] = {(i0+i1)/2.0f, (j0+j1)/2.0f};
  float S_m[4] = {
    ((i0-mu_m[0])*(i0-mu_m[0]) + (i1-mu_m[0])*(i1-mu_m[0]))/2.0f,
    ((i0-mu_m[0])*(j0-mu_m[1]) + (i1-mu_m[0])*(j1-mu_m[1]))/2.0f,
    ((i0-mu_m[0])*(j0-mu_m[1]) + (i1-mu_m[0])*(j1-mu_m[1]))/2.0f,
    ((j0-mu_m[1])*(j0-mu_m[1]) + (j1-mu_m[1])*(j1-mu_m[1]))/2.0f};

  // again convolve:
  for(int k=0;k<2;k++) mu[k] += mu_m[k];
  for(int k=0;k<4;k++) S [k] += S_m [k];
#endif

  // TODO: depth of field?
#if 1
  // approximate circle of confusion?
#endif

  // do not use offset:
  for(int k=0;k<2;k++) mu[k] = 0.0f;

  float col[3] = {0.0f};
  for(int i=0;i<num_wavelengths;i++)
  {
    // perpath.lambda = lam[i];
    if(!(val[i] > 0.0) || !(pdf[i] > 0.0)) continue;
    double v = val[i] / pdfsum  * mis[i];
    if(!(v == v)) continue;
    float c2[3];
    spectrum_p_to_camera(lam[i], v, c2);
    for(int k=0;k<3;k++) col[k] += c2[k];
  }
  // view_splat_col(&perpath, col);
  view_splat_gaussian(&perpath, mu, S, col);
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

void pointsampler_prepare_frame(pointsampler_t *s)
{
  if(rt.threads->end >> 32 > s->reinit)
    halton_init_random(&s->h, rt.anim_frame + ++s->reinit);
}

void pointsampler_stop_learning(pointsampler_t *s) {}
