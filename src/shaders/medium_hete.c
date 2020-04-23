/*
    This file is part of corona-13.

    copyright (c) 2015--2018 johannes hanika.

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

// heterogeneous medium, wrapping around custom voxel format.

#include "vol/vol.h"
#include "vol/trace.h"
#include "vol/interpolation.h"
#include "vol/lighthierarchy.h"
#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum.h"
#include "prims.h"
#include "display.h"
#include "lights.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FNEE_LOD 1
#define FNEE_USE_MAX 0

typedef struct
{
  // mean cosine
  float g0, g1;
  float sigma_e; // to get mu_e = density * sigma_e
  float sigma_t; // to get mu_t = density * sigma_t
  float sigma_s; // to get mu_s = density * sigma_s
  vol_interpolation_t interpolation;

  vol_tree_t *tree;
  int mshader;
}
medium_t;


// andrea says these are the right params for clouds:
// > R:{0.943289,-0.272228,0.821596}
// > 
// > G: {0.952029,-0.289216,0.812471}
// > 
// > B: {0.956937,-0.342329,0.81399} 
// > 
// > Values are for RGB g1, g2, and weight, i.e. pHG(Cos[t], fit[0])*fit[2]
// > + pHG(Cos[t], fit[1])*(1-fit[2]) 
// > 
// > If you are after something specific like non-mixed, larger or smaller,
// > let me know.

static inline void mie_fit(
    medium_t *s,
    const float lambda,
    float *g1,
    float *g2,
    float *t)
{
  *g1 = s->g0;
  *g2 = s->g1;
  *t = .5f;
#if 0
  const float fit[3][3] = {
    {0.943289,-0.272228,0.821596},
    {0.952029,-0.289216,0.812471},
    {0.956937,-0.342329,0.81399}};
  const float l[3] = {400.0, 550.0, 650.0};
  float tt = 0.0;
  int w = 0;
  if(lambda > l[0] && lambda <= l[1])
  {
    w = 0;
    tt = (lambda - l[0])/(l[1]-l[0]);
  }
  else if(lambda > l[1] && lambda <= l[2])
  {
    w = 1;
    tt = (lambda - l[1])/(l[2]-l[1]);
  }
  else if(lambda > l[2])
  {
    w = 1;
    tt = 1;
  }

  *g1 = (1.0f-tt)*fit[w][0] + tt*fit[w+1][0];
  *g2 = (1.0f-tt)*fit[w][1] + tt*fit[w+1][1];
  *t  = (1.0f-tt)*fit[w][2] + tt*fit[w+1][2];
#endif
}

static inline float mean_cos(medium_t *s)
{
  return (s->g0 + s->g1)*.5f; // blatant lie
  // return 0.9; // approximately right for the mie fit above, maybe
}

// only use interpolated density on coarser levels (faster)
// #define INTERP(s, e) (_lod(e) ? s->interpolation : (s->interpolation & ~s_vol_trilinear))
// #define INTERP(s, e) (s->interpolation & ~s_vol_trilinear)
#define INTERP(s, e) s->interpolation

static inline float regularise(int plen)
{
  return 1.0; // no similarity relation
  float start = 2, end = 3;
  const float t = CLAMP((plen - start)/(end-start), 0.0f, 1.0f);
  return 1.-t;
}

// switch on volume callbacks, to distinguish between hete and homo volumes:
int volume_enabled(void *data)
{
  return 1;
}

mf_t volume_transmittance(path_t *p, int e, void *data)
{
  medium_t *s = data;
  mf_t emission = mf_set1(0.0f);
  const float reg = regularise(p->length);
  const float g = mean_cos(s);
  const float sigma_reg = (1-g)*(1-reg)+reg;
#ifdef TRANS_TRACK_LEN
  const mf_t transmittance = vol_trace_transmittance_track_length(
#elif defined(TRANS_RATIO_TRACK)
  const mf_t transmittance = vol_trace_transmittance_ratio_tracking(
#elif defined(TRANS_RES_RATIO)
  const mf_t transmittance = vol_trace_transmittance_residual_ratio_tracking(
#else
  const mf_t transmittance = vol_trace_transmittance(
#endif
      s->tree, p->v[e-1].hit.x, p->e[e].omega,
      p->e[e].dist, mf_set1(sigma_reg*s->sigma_t), INTERP(s, e), 0, p->time,
      p->lambda, s->sigma_e > 0.0f ? &emission : 0, 0);
#ifdef SEGMENT_EMISSION
  if(s->sigma_e > 0.0f) p->e[e].contribution = mf_mul(mf_set1(sigma_reg*s->sigma_e), emission);
#endif
  p->e[e].transmittance = transmittance;
  p->e[e].pdf = mf_set1(1.0f); // was deterministic
  return transmittance;
}

float volume_sample(path_t *p, int e, void *data)
{
  medium_t *s = data;
  p->e[e].contribution = mf_set1(0.0f);
  mf_t emission = mf_set1(0.0f), transmittance = mf_set1(1.0f);
  float last_density = -1.0f;
  const float reg = regularise(p->length);
  const float g = mean_cos(s);
  const float sigma_reg = (1-g)*(1-reg)+reg;
  const float dist = vol_trace_sample(s->tree, p->v[e-1].hit.x, p->e[e].omega,
      p->e[e].dist, mf_set1(sigma_reg*s->sigma_t), pointsampler(p, s_dim_free_path), INTERP(s, e), 0, p->time,
      p->lambda, s->sigma_e > 0.0f ? &emission : 0, &transmittance, 0, &last_density);
#ifdef SEGMENT_EMISSION
  if(s->sigma_e > 0.0f) p->e[e].contribution = mf_mul(mf_set1(sigma_reg * s->sigma_e), emission);
#endif
  p->e[e].pdf = p->e[e].transmittance = transmittance;
#if 0
  const float mu_t = sigma_reg * s->sigma_t * last_density; // TODO: replace by something numerically stable and fast!
#else // :( need that for numeric precision on voxel boundaries
  float dens_temp[2] = {0.0f};
  vol_lookup(s->tree,
      p->v[e-1].hit.x[0] + MIN(p->e[e].dist, dist) * p->e[e].omega[0],
      p->v[e-1].hit.x[1] + MIN(p->e[e].dist, dist) * p->e[e].omega[1],
      p->v[e-1].hit.x[2] + MIN(p->e[e].dist, dist) * p->e[e].omega[2],
      s_vol_density | s_vol_temperature, INTERP(s, v), 0, p->time, dens_temp);
  const mf_t mu_t = mf_set1(sigma_reg * s->sigma_t * dens_temp[0]);
#endif
  if(dist < p->e[e].dist)
    p->e[e].pdf = mf_mul(p->e[e].pdf, mu_t);
  // if(dist <= p->e[e].dist)
    p->v[e].interior.mu_t = mu_t;
  // else
    // p->v[e].interior.mu_t = 0;
  return dist;
}

mf_t volume_pdf_adj(const path_t *p, int e, void *data)
{
  medium_t *s = data;
  const float reg = regularise(p->length);
  const float g = mean_cos(s);
  const float sigma_reg = (1-g)*(1-reg)+reg;
  if(mf(p->e[e].transmittance,0) > 0.0f) // assume good cached value
  {
    if(p->v[e-1].material_modes & s_volume) return mf_mul(mf_mul(p->e[e].transmittance, mf_set1(sigma_reg)), p->v[e-1].interior.mu_t);
    return p->e[e].transmittance;
  }
  float last_density = 0.0f;
  float dist = p->e[e].dist;
  if(p->v[e].flags & s_environment)
  { // compute real distance between the points:
    float delta[3] = {
      p->v[e].hit.x[0] - p->v[e-1].hit.x[0],
      p->v[e].hit.x[1] - p->v[e-1].hit.x[1],
      p->v[e].hit.x[2] - p->v[e-1].hit.x[2]};
    dist = sqrtf(dotproduct(delta, delta));
  }
  float wi[3] = { -p->e[e].omega[0], -p->e[e].omega[1], -p->e[e].omega[2] };
  const mf_t transmittance = vol_trace_transmittance(s->tree, p->v[e].hit.x, wi,
      dist, mf_set1(sigma_reg * s->sigma_t), INTERP(s, e), 0, p->time,
      p->lambda, 0, &last_density);
  if(p->v[e-1].material_modes & s_volume) return mf_mul(transmittance, mf_mul(mf_set1(sigma_reg), p->v[e-1].interior.mu_t));
  return transmittance;
}

mf_t volume_pdf(const path_t *p, int e, void *data)
{
  medium_t *s = data;
  const float reg = regularise(p->length);
  const float g = mean_cos(s);
  const float sigma_reg = (1-g)*(1-reg)+reg;
  if(mf(p->e[e].transmittance,0) > 0.0f) // assume good cached value
  {
    if(p->v[e].material_modes & s_volume) return
      mf_mul(mf_mul(p->e[e].transmittance, mf_set1(sigma_reg)), p->v[e].interior.mu_t);
    return p->e[e].transmittance;
  }
  const mf_t transmittance = vol_trace_transmittance(s->tree, p->v[e-1].hit.x, p->e[e].omega,
      p->e[e].dist, mf_set1(sigma_reg * s->sigma_t), INTERP(s, e), 0, p->time,
      p->lambda, 0, 0);
  if(p->v[e].material_modes & s_volume)
    return mf_mul(transmittance, mf_mul(mf_set1(sigma_reg), p->v[e].interior.mu_t));
  return transmittance;
}

// called after init, in case we are somehow associated with a geo object.
// (do this by assigning the shape our shaderid)
// we interpret this to rescale to the object's bounding box.
int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = self->data;
  s->mshader = self - rt.shader->shader;
  if(s->mshader != rt.shader->exterior_medium_shader) return 0; // we want to rescale and discard only if we're the exterior medium
  float aabb[6];
  prims_get_shape_aabb(rt.prims, shapeid, aabb);
  fprintf(stderr, "[medium_hete] captured shape[%d] (%s) fitting vol %gx%gx%g into bounding box %gx%gx%g\n",
      shapeid, rt.prims->shape[shapeid].name,
      s->tree->content_box[3]-s->tree->content_box[0],
      s->tree->content_box[4]-s->tree->content_box[1],
      s->tree->content_box[5]-s->tree->content_box[2],
      aabb[3]-aabb[0], aabb[4]-aabb[1], aabb[5]-aabb[2]);

  // we want to count against scene bounds
  prims_extend_ghost_aabb(rt.prims, aabb);

  // TODO: allow rotations?
  // adjust voxel size to new scale, only use first dimension (rest isotropic scaling)
  const float scale = (aabb[3]-aabb[0])/(s->tree->content_box[3]-s->tree->content_box[0]);
  s->tree->voxel_size = s->tree->header->voxel_size * scale;
  float loc[3];
  for(int k=0;k<3;k++) loc[k] = aabb[k]-s->tree->content_box[k]*scale;
  for(int k=0;k<6;k++) s->tree->content_box[k] *= scale;
  for(int k=0;k<6;k++) s->tree->aabb[k] *= scale;
  vol_create_transform(s->tree, loc, 0);

  return 1; // we want to discard the shape (only used for proxy aabb)
}

float prepare(path_t *p, int v, void *data)
{ 
  medium_t *s = data;
  // set volume properties on vertex
  p->v[v].interior.mean_cos = regularise(p->length) * mean_cos(s);
  p->v[v].interior.shader = s->mshader;
  p->v[v].shading.em = mf_set1(0.0f);
  // envmap vertices are never prepare()d, so invalid primid means scattering:
  if(primid_invalid(p->v[v].hit.prim))
  {
    const float reg = regularise(p->length);
    const float g = mean_cos(s);
    const float sigma_reg = (1-g)*(1-reg)+reg;
    // query density and derive mu_t, mu_s from that.
    float dens_temp[2] = {0.0f};
    vol_lookup(s->tree, p->v[v].hit.x[0], p->v[v].hit.x[1], p->v[v].hit.x[2],
        s_vol_density | s_vol_temperature, INTERP(s, v), 0, p->time, dens_temp);
    p->v[v].interior.mu_t = mf_set1(sigma_reg * dens_temp[0] * s->sigma_t);
    p->v[v].interior.mu_s = mf_set1(sigma_reg * dens_temp[0] * s->sigma_s);
    
    // set material modes only if we are surely inside the volume
    if (dens_temp[0] > 0)
      p->v[v].material_modes = s_volume | s_glossy;

#ifndef SEGMENT_EMISSION // can't use point emission and line emission at the same time:
    if(s->sigma_e > 0.0f)
    {
      if(dens_temp[1] > 0.0f && dens_temp[0] > 0.0f)
      {
        p->v[v].material_modes = p->v[v].mode = p->v[v].mode | s_volume | s_emit | s_diffuse;
        // L_e = mu_e * L [W/m^3] = sigma_e * density * shader
        p->v[v].shading.em = mf_mul(mf_set1(sigma_reg * dens_temp[0] * s->sigma_e), s->tree->shader(dens_temp[0], dens_temp[1], p->lambda));
      }
    }
#else
    p->v[v].shading.em = mf_set1(0.0f);
    if(s->sigma_e > 0.0f)
      if(dens_temp[1] > 0.0f && dens_temp[0] > 0.f)
        p->v[v].material_modes =
          p->v[v].mode = p->v[v].mode | s_volume | s_emit | s_diffuse;
#endif
  }
  return 1.0f;
}

mf_t sample(path_t *p, void *data)
{
  medium_t *s = data;
  const int v = p->length-1;
  if(mf_all(mf_lte(p->v[v].interior.mu_t, mf_set1(0.0f)))) return mf_set1(0.0f); // kill paths in the case that we sampled a position in vacuum
  p->v[v].mode |= s_glossy | s_volume;
  const hit_t *hit = &p->v[v].hit;
  float out[3];
  float g1, g2, t;
  mie_fit(s, mf(p->lambda,0), &g1, &g2, &t);
  const float reg = regularise(p->length);
  float g = g2 * reg;
  if(pointsampler(p, s_dim_scatter_mode) < t) g = g1 * reg;
  sample_hg(g, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), out, 0);
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = hit->n[k] * out[0] + hit->a[k] * out[1] + hit->b[k] * out[2];
  p->v[v+1].pdf = mf_set1(
       t * sample_eval_hg(g1 * reg, p->e[v].omega, p->e[v+1].omega) +
    (1-t)* sample_eval_hg(g2 * reg, p->e[v].omega, p->e[v+1].omega));
  g = mean_cos(s);
  return mf_mul(mf_set1((1-g)*(1-reg)+reg), p->v[v].interior.mu_s);
}

mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  medium_t *s = data;
  if(!(p->v[v].mode & s_volume)) return mf_set1(0.0f);
  const float reg = regularise(p->length);
  float g1, g2, t;
  mie_fit(s, mf(p->lambda,0), &g1, &g2, &t);
  return mf_set1(
       t * sample_eval_hg(g1 * reg, p->e[e1].omega, p->e[e2].omega) +
    (1-t)* sample_eval_hg(g2 * reg, p->e[e1].omega, p->e[e2].omega));
}

mf_t brdf(path_t *p, int v, void *data)
{
  medium_t *s = data;
  p->v[v].mode |= s_glossy | s_volume;
  if(mf_all(mf_lte(p->v[v].interior.mu_t, mf_set1(0.0f))))
    return mf_set1(0.0f); // kill paths in the case that we sampled a position in vacuum
  const float reg = regularise(p->length);
  const float g = mean_cos(s);
  float g1, g2, t;
  mie_fit(s, mf(p->lambda,0), &g1, &g2, &t);
  return mf_mul(p->v[v].interior.mu_s, 
    mf_set1(((1-g)*(1-reg)+reg) *
    (   t *sample_eval_hg(g1 * reg, p->e[v].omega, p->e[v+1].omega) +
     (1-t)*sample_eval_hg(g2 * reg, p->e[v].omega, p->e[v+1].omega))));
}

int init(FILE* f, void** data)
{
  medium_t *s = malloc(sizeof(*s));
  memset(s, 0, sizeof(medium_t));
  s->interpolation = s_vol_constant;
  // s->interpolation = s_vol_smooth;
  *data = (void *)s;
  char filename[1024];
  s->sigma_t = 10.0f;
  s->sigma_s = 10.0f;
  s->sigma_e = 0.0;
  int i = fscanf(f, "%f %f %f %f %f %s", &s->g0, &s->g1, &s->sigma_s, &s->sigma_t, &s->sigma_e, filename);
  if(i != 6)
  {
    fprintf(stderr, "[medium_hete] could not parse all arguments! expecting medium_hete <mean_cosine1> <mean_cosine2> <sigma_s> <sigma_t> <sigma_e> <tree_filename>\n");
    s->g0 = 0.0f;
    s->g1 = 0.0f;
    return 1;
  }
  int res = fscanf(f, "%*[^\n]\n");

  display_control_add(rt.display, "[hete] extc", &s->sigma_t, 0.0001, 1000.0, 2.0, 1, 1);
  display_control_add(rt.display, "[hete] scat", &s->sigma_s, 0.0001, 1000.0, 2.0, 1, 1);
  display_control_add(rt.display, "[hete] emit", &s->sigma_e, 0.0, 1000.0, 1.0, 0, 1);
  display_control_add(rt.display, "[hete] g0", &s->g0, -1, 1, 0.1, 0, 1);
  display_control_add(rt.display, "[hete] g1", &s->g1, -1, 1, 0.1, 0, 1);

  s->tree = vol_open(filename);
  if(!s->tree)
  {
    char fname[1024];
    snprintf(fname, 1024, "%s/%s", rt.searchpath, filename);
    s->tree = vol_open(fname);
  }
  if(!s->tree)
    fprintf(stderr, "[medium_hete] could not open volume data `%s'!\n", filename);

  fprintf(stderr, "[medium_hete] opened volume with %luM voxels\n", vol_count_voxels(s->tree)/1000000);

  return res == -1;
}

void cleanup(void *data)
{
  medium_t *s = data;
  vol_close(s->tree);
}

mf_t volume_pdf_nee(const path_t *p, int v, void *data)
{
  medium_t *s = data;
  // XXX FIXME: repair this to be hwl aware!
  return vol_lighthierarchy_pdf_point(
      s->tree, p->lambda, p->time, 0, 0,
      (p->v[v-1].material_modes & (s_volume | s_fiber)) ?
      0 : p->v[v-1].hit.x,
      p->v[v-1].hit.n,
      p->v[v].hit.x);
}

/* sample end point v for next event estimation */
mf_t volume_sample_nee(path_t *p, int v, void *data)
{
  medium_t *s = data;
  const mf_t L_e = mf_mul(mf_set1(s->sigma_e), vol_lighthierarchy_sample_point(
      s->tree, p->lambda, p->time,
      0,
      0,
      pointsampler(p, s_dim_nee_x),
      pointsampler(p, s_dim_nee_y),
      pointsampler(p, s_dim_nee_light2),
      (p->v[v-1].material_modes & (s_volume | s_fiber)) ?
      0 : p->v[v-1].hit.x,
      p->v[v-1].hit.n,
      p->v[v].hit.x,
      &p->v[v].pdf)); // is in volume area measure 1/m^3

#if 0 // needed for lod i think:
  float dens_temp[2] = {0.0f};
  vol_lookup(s->tree, p->v[v].hit.x[0], p->v[v].hit.x[1], p->v[v].hit.x[2],
      s_vol_density | s_vol_temperature, INTERP(s, v), 0, p->time, dens_temp);
  const float L_e = s->sigma_e * dens_temp[0] * s->tree->shader(dens_temp[0], dens_temp[1], p->lambda);
#endif

  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].mode = p->v[v].material_modes = s_volume | s_diffuse | s_emit;
  p->v[v].shading.em = L_e;
  p->v[v].material_modes = s_volume | s_diffuse | s_emit;
  if(mf(p->v[v].pdf,0) <= 0.0f) return mf_set1(0.0f);
  return mf_div(L_e, p->v[v].pdf);
}

mf_t volume_pdf_fnee_direction(const path_t *p, int e, void *data)
{
  medium_t *s = data;
  return vol_lighthierarchy_pdf_direction(
      s->tree,
      p->v[e-1].hit.x,
      p->e[e].omega,
      p->lambda, p->time,
      s->interpolation,
      FNEE_LOD,
      FNEE_USE_MAX,
      (p->v[e-1].material_modes & (s_volume | s_fiber)) ?
      0 : p->v[e-1].hit.x,
      p->v[e-1].hit.n);
}

/* sample direction p->e[e].{omega,pdf} for forward next event estimation */
void volume_sample_fnee_direction(path_t *p, int e, void *data)
{
  medium_t *s = data;
  float pos[3] = {0.0f};

  vol_lighthierarchy_sample_point(
      s->tree, p->lambda, p->time,
      FNEE_LOD,
      FNEE_USE_MAX,
      pointsampler(p, s_dim_nee_x),
      pointsampler(p, s_dim_nee_y),
      pointsampler(p, s_dim_nee_light2),
      (p->v[e-1].material_modes & (s_volume | s_fiber)) ?
      0 : p->v[e-1].hit.x,
      p->v[e-1].hit.n,
      pos,
      0);

  for(int k=0;k<3;k++)
    p->e[e].omega[k] = pos[k] - p->v[e-1].hit.x[k]; 
  p->e[e].dist = sqrtf(dotproduct(p->e[e].omega, p->e[e].omega));
  for(int k=0;k<3;k++)
    p->e[e].omega[k] /= p->e[e].dist;

  p->v[e].pdf = volume_pdf_fnee_direction(p, e, data);
}

mf_t volume_pdf_fnee_direction_adj(const path_t *p, int e, void *data)
{
  return volume_pdf_fnee_direction(p, e, data);
}

/* sample start point v[0] for light tracing */
mf_t volume_sample_light(path_t *p, void *data)
{
	assert(0 && "not implemented");
	return mf_set1(0.f);
}

mf_t volume_pdf_light(const path_t *p, int v, void *data)
{
	assert(0 && "not implemented");
	return mf_set1(0.f);
}
