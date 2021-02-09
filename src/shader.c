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
#include "prims.h"
#include "pathspace/manifold.h"
#include "render.h"
#include "shader.h"
#include "sampler_common.h"
#include "sampler.h"
#include "spectrum.h"
#include "lights.h"
#include "accel.h"
#include "shaders/daylight.h"

#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int shader_vol_hete(const path_t *p, int e)
{
  if(p->e[e].vol.shader < 0) return 0; // no volume at all
  shader_so_t *s = rt.shader->shader + p->e[e].vol.shader;
  if(s->volume_enabled && s->volume_enabled(s->data))
    return 1;
  return 0; // homogeneous volume
}

// volume related things on the edge.
// compute transmittance, potentially write edge emission
mf_t shader_vol_transmittance(path_t *p, int e)
{
  p->e[e].contribution = mf_set1(0.0f);
  p->e[e].pdf = mf_set1(1.0f);
  if(p->e[e].vol.shader >= 0)
  {
    shader_so_t *s = rt.shader->shader + p->e[e].vol.shader;
    if(s->volume_transmittance)
    {
      p->e[e].transmittance = s->volume_transmittance(p, e, s->data);
    }
    else
    { // default case: homogeneous scattering.
      // clamp environment distance at 1km so sun will still have an influence..
      if((p->v[e].flags & s_environment) || (p->v[e-1].flags & s_environment))
        p->e[e].transmittance = mf_exp(mf_mul(mf_set1(-10000.0f), p->e[e].vol.mu_t));
      else
        p->e[e].transmittance = mf_exp(mf_mul(mf_set1(-p->e[e].dist), p->e[e].vol.mu_t));
      p->e[e].pdf = mf_set1(1.0f);
      p->e[e].contribution = mf_set1(0.0f);
      return p->e[e].transmittance;
    }
  }
  else p->e[e].transmittance = mf_set1(1.0f);
  return p->e[e].transmittance;
}

// sample free distance
float shader_vol_sample(path_t *p, int e)
{
  float dist = FLT_MAX;
  p->e[e].contribution = mf_set1(0.0f);
  p->e[e].pdf = mf_set1(1.0f);
  p->e[e].transmittance = mf_set1(1.0f);

  if(p->e[e].vol.shader >= 0)
  { // in medium
    shader_so_t *s = rt.shader->shader + p->e[e].vol.shader;
    if(s->volume_sample)
    {
      return s->volume_sample(p, e, s->data);
    }
    else
    { // default case: homogeneous scattering.
      if(mf(p->e[e].vol.mu_s, 0) > .0f)
      { // in scattering medium
        const float rf = pointsampler(p, s_dim_free_path);
        dist = - logf(1.0f - rf)/mf(p->e[e].vol.mu_t, 0);
        if(!(dist > 0.0)) dist = 1e-15;
        p->e[e].pdf = p->e[e].transmittance = mf_exp(mf_mul(mf_set1(-dist), p->e[e].vol.mu_t));
        if(dist < p->e[e].dist) p->e[e].pdf = mf_mul(p->e[e].pdf, p->e[e].vol.mu_t);
      }
      else
        p->e[e].transmittance = mf_exp(mf_mul(mf_set1(-p->e[e].dist), p->e[e].vol.mu_t));
    }
  }
  assert(dist > 0.0);
  return dist;
}

// pdf for free path sampling method.
mf_t shader_vol_pdf(const path_t *p, int e)
{
  mf_t pdf = mf_set1(1.0f);
  if(p->e[e].vol.shader >= 0)
  {
    shader_so_t *s = rt.shader->shader + p->e[e].vol.shader;
    if(s->volume_pdf)
    {
      return s->volume_pdf(p, e, s->data);
    }
    else
    { // default case: homogeneous scattering.
      // pdf only changes in scattering medium
      if(mf(p->e[e].vol.mu_s, 0) > .0f)
      {
        pdf = mf_exp(mf_mul(mf_set1(-p->e[e].dist), p->e[e].vol.mu_t));
        // need to mul mu_t only if v[v] is a medium interaction
        if(!(p->v[e].flags & s_environment) &&
           primid_invalid(p->v[e].hit.prim)) pdf = mf_mul(pdf, p->e[e].vol.mu_t);
      }
    }
  }
  return pdf;
}

mf_t shader_vol_pdf_adjoint(const path_t *p, int e)
{
  mf_t pdf = mf_set1(1.0f);
  if(p->e[e].vol.shader >= 0)
  { // in medium
    shader_so_t *s = rt.shader->shader + p->e[e].vol.shader;
    if(s->volume_pdf_adj)
    {
      return s->volume_pdf_adj(p, e, s->data);
    }
    else
    { // homogeneous
      if(mf(p->e[e].vol.mu_s, 0) > .0f)
      {
        pdf = mf_exp(mf_mul(mf_set1(-p->e[e].dist), p->e[e].vol.mu_t));
        if(!(p->v[e-1].flags & s_environment) &&
           primid_invalid(p->v[e-1].hit.prim)) pdf = mf_mul(pdf, p->e[e].vol.mu_t);
      }
    }
  }
  return pdf;
}

float prepare_d(path_t *p, int v, void *data)
{
  if(mf_any(mf_gt(p->v[v].shading.rd, mf_set1(0.0f))))
    p->v[v].material_modes = s_reflect | s_diffuse;
  return 1.0f;
}

// default diffuse white shader: brdf = rd/pi, p=cos/pi
mf_t sample_d(path_t *p, void *data)
{
  const int v = p->length; // about to sample that vertex

  const float x1 = pointsampler(p, s_dim_omega_x);
  const float x2 = pointsampler(p, s_dim_omega_y);
  float s = sqrtf(x1);
  // light tracer samples geometric normal
  float *n = (p->v[0].mode & s_emit) ? p->v[v-1].hit.gn : p->v[v-1].hit.n;
  float sign = 1.0f;
  if((p->v[0].mode & s_emit) && (p->v[v-1].flags & s_inside))
    sign = -1.0f; // need to flip geo normal, too.
  for(int k=0;k<3;k++)
    p->e[v].omega[k] =
      sqrtf(1.0 - x1)   * sign * n[k] +
      s*cosf(2*M_PI*x2) * p->v[v-1].hit.a[k] +
      s*sinf(2*M_PI*x2) * p->v[v-1].hit.b[k];
  p->v[v].pdf = mf_set1(1.0f/M_PI);

  const float cos_out_ng = dotproduct(p->v[v-1].hit.gn, p->e[v].omega);
  if(p->v[v-1].flags & s_inside)
  {
    if(cos_out_ng >= 0.0f) return mf_set1(0.0f);
  }
  // XXX
  // else if(cos_out_ng <= 0.0f) return mf_set1(0.0f);

  mf_t throughput = mf_set1(0.0f);
  if(p->v[0].mode & s_emit)
  {
    const float cos_ns = dotproduct(p->v[v-1].hit.n, p->e[v-1].omega);
    const float cos_ng = dotproduct(p->v[v-1].hit.gn, p->e[v-1].omega);
    // undo sampling the geometric cosine during light tracing, clamp to avoid excessive variance:
    throughput = mf_mul(p->v[v-1].shading.rd, mf_set1(fminf(4.0f, fabsf(cos_ns/cos_ng))));
  }
  else throughput = p->v[v-1].shading.rd;

  // only set mode after we're sure we will return > 0:
  if(mf_any(mf_gt(throughput, mf_set1(0.0f)))) p->v[v-1].mode = s_diffuse | s_reflect;
  return throughput;
}

mf_t brdf_d(path_t *p, int v, void *data)
{
  p->v[v].mode = s_diffuse | s_reflect;
  const float cos_out_ns = dotproduct(p->v[v].hit.n, p->e[v+1].omega);
  if(cos_out_ns <= 0) return mf_set1(0.0f);
  return mf_mul(p->v[v].shading.rd, mf_set1(1.0f/M_PI)); // XXX
  const float cos_out_ng = dotproduct(p->v[v].hit.gn, p->e[v+1].omega);
  const float cos_in_ns = -dotproduct(p->v[v].hit.n, p->e[v].omega);
  if(p->v[0].mode & s_emit)
  {
    // shading normal madness. we need two ratios since path space uses
    // the shading normal hit.n in the geometric term, not the geometric one.
    const float cos_in_ng = dotproduct(p->v[v].hit.gn, p->e[v].omega);
    // if(cos_in_ns > 0.0f && cos_out_ns > 0.0f) // reciprocal but black borders.
    if((!(p->v[v].flags & s_inside) && (cos_out_ng > 0.0f)) || ((p->v[v].flags & s_inside) && (cos_out_ng < 0.0f)))
      return mf_mul(mf_set1(fminf(4.0f, fabsf(cos_in_ns*cos_out_ng/(cos_in_ng * cos_out_ns)))/M_PI), p->v[v].shading.rd);
  }
  // else if(cos_in_ns > 0.0f && cos_out_ns > 0.0f) // reciprocal but black borders.
  else if(cos_out_ns > 0.0f) // ignore incoming under the surface
  // else
  {
    if(p->v[v].flags & s_inside)
    {
      if(cos_out_ng >= 0.0f) return mf_set1(0.0f);
    }
    else if(cos_out_ng <= 0.0f) return mf_set1(0.0f);
    return mf_mul(p->v[v].shading.rd, mf_set1(1.0f/M_PI));
  }
  return mf_set1(0.0f);
}

mf_t pdf_d(const path_t *p, int e1, int v, int e2, void *data)
{
  return mf_set1(1.0f/M_PI);
}

extern int init_d(FILE* f, void** data)
{
  *data = NULL;
  int dreggn = fscanf(f, "%*[^\n]\n");
  if(dreggn == -1) fprintf(stderr, "gcc stinks\n");
  return 0;
}

mf_t sky_black(const path_t *p, int v, void* data)
{
  return mf_set1(0.0f);
}

mf_t sky_cloudy(const path_t *p, int v, void* data)
{
  const float power = 1.0f;
  const float scale = 500.0f;

  if(v == 0)
    return mf_set1(power * scale * 0.5f*(1.0 - p->e[1].omega[2]));
  else
    return mf_set1(power * scale * 0.5f*(1.0 + p->e[v].omega[2]));
}

mf_t sky_cloudy_sample(path_t *p, void *data)
{
  const float power = 1.0f;
  const float scale = 500.0f;

  float x1, x2;
  if(p->v[0].mode & s_emit)
  { // light ray
    x1 = pointsampler(p, s_dim_edf_x);
    x2 = pointsampler(p, s_dim_edf_y);
  }
  else
  { // next event
    x1 = pointsampler(p, s_dim_nee_x);
    x2 = pointsampler(p, s_dim_nee_y);
  }
  // sample cos on full sphere, pdf = (.5+z/2) * 1/(2 pi)
  const float z = -(1.0f - 2.0f*sqrtf(1.0f - x1));
  const float sin_theta = sqrtf(1.0 - z*z);
  const float x = sin_theta * cosf(2.f*M_PI*x2);
  const float y = sin_theta * sinf(2.f*M_PI*x2);

  const int v = p->length; // sample new vertex at the end
  p->v[v].shading.em = mf_set1((.5f+z*.5f)*power * scale);
  p->v[v].pdf = mf_set1((.5f + z*.5f)/(2.0f * M_PI));
  p->v[v].flags = s_environment;
  p->v[v].mode = s_emit;
  p->v[v].shading.roughness = 1.0f;
  p->e[v].omega[2] = z;
  p->e[v].omega[0] = x;
  p->e[v].omega[1] = y;
  p->e[v].dist = FLT_MAX;
  if(p->length)
  { // only set this in case we're doing next event estimation
    const float *aabb = accel_aabb(rt.accel);
    const float far = aabb[3] + aabb[4] + aabb[5]
                    - aabb[0] - aabb[1] - aabb[2];
    for(int k=0;k<3;k++)
    {
      p->v[v].hit.x[k] = p->v[v-1].hit.x[k] + far*p->e[v].omega[k];
      p->v[v].hit.n[k] = p->v[v].hit.gn[k] = - p->e[v].omega[k];
    }
  }
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = -1;
  return mf_div(p->v[v].shading.em, p->v[v].pdf);
}

mf_t sky_cloudy_pdf(const path_t *p, int v, void *data)
{
  return mf_set1((0.5f + p->e[v].omega[2]*.5f)/(2.0f * M_PI));
}

mf_t shader_sky_sample_d(path_t *p, void *data)
{
  const int v = p->length;
  p->v[v].mode = s_absorb;
  p->v[v].throughput = mf_set1(0.0f);
  p->v[v].shading.em = mf_set1(0.0f);
  return mf_set1(0.0f);
}

mf_t shader_sky_pdf_d(const path_t *p, int v, void *data)
{
  return mf_set1(0.0f);
}

mf_t shader_sky_sample(path_t *p)
{
  const mf_t throughput = rt.shader->skyshader.sample(p, rt.shader->skyshader.data);
  // also sample the point position out of the bounding box
  const float x1 = pointsampler(p, s_dim_light_x);
  const float x2 = pointsampler(p, s_dim_light_y);

  // abused next event estimation code, make it a light ray now.
  for(int k=0;k<3;k++)
  {
    p->v[0].hit.gn[k] =
    p->v[0].hit.n[k]  =
    p->e[1].omega[k]  = -p->e[0].omega[k];
  }

  get_onb(p->v[0].hit.gn, p->v[0].hit.a, p->v[0].hit.b);

  // create quad large enough to cover aabb
  const float *aabb = accel_aabb(rt.accel);
  float mina = INFINITY, maxa = - INFINITY, minb = INFINITY, maxb = - INFINITY;
  for(int i=0;i<4;i+=3) for(int j=0;j<4;j+=3) for(int k=0;k<4;k+=3) 
  {
    const float dota = aabb[i+0]*p->v[0].hit.a[0] + aabb[j+1]*p->v[0].hit.a[1] + aabb[k+2]*p->v[0].hit.a[2];
    const float dotb = aabb[i+0]*p->v[0].hit.b[0] + aabb[j+1]*p->v[0].hit.b[1] + aabb[k+2]*p->v[0].hit.b[2];
    if(dota < mina) mina = dota;
    if(dota > maxa) maxa = dota;
    if(dotb < minb) minb = dotb;
    if(dotb > maxb) maxb = dotb;
  }

  // for p(pos) = 1/A
  const float pdf_area = 1.0f/((maxa-mina)*(maxb-minb));
  p->v[1].pdf = mf_mul(p->v[0].pdf, mf_set1(pdf_area)); // vertex is sampled in vertex area measure on the quad
  p->v[0].pdf = mf_set1(1.0f);
    
  const float far = dotproduct(p->v[0].hit.gn, aabb)
                    + aabb[3] + aabb[4] + aabb[5]
                    - aabb[0] - aabb[1] - aabb[2];

  for(int k=0;k<3;k++)
    p->v[0].hit.x[k] = - far*p->v[0].hit.gn[k] +
      p->v[0].hit.a[k]*(mina + (maxa - mina)*x1) +
      p->v[0].hit.b[k]*(minb + (maxb - minb)*x2);

  return mf_div(throughput, mf_set1(pdf_area));
}

mf_t shader_sky_pdf(const path_t *p, int v)
{
  mf_t pdf = rt.shader->skyshader.pdf(p, v, rt.shader->skyshader.data);

  assert(v > 0); // in that case please call next event instead
  float a[3], b[3];
  // we are called with v=1 for the adjoint case and v=length-1 for the same direction.
  // in both cases v is the edge leading up to vertex v
  get_onb(p->e[v].omega, a, b);

  // create quad large enough to cover aabb
  const float *aabb = accel_aabb(rt.accel);
  float mina = INFINITY, maxa = - INFINITY, minb = INFINITY, maxb = - INFINITY;
  for(int i=0;i<4;i+=3) for(int j=0;j<4;j+=3) for(int k=0;k<4;k+=3) 
  {
    const float dota = aabb[i+0]*a[0] + aabb[j+1]*a[1] + aabb[k+2]*a[2];
    const float dotb = aabb[i+0]*b[0] + aabb[j+1]*b[1] + aabb[k+2]*b[2];
    if(dota < mina) mina = dota;
    if(dota > maxa) maxa = dota;
    if(dotb < minb) minb = dotb;
    if(dotb > maxb) maxb = dotb;
  }

  // for p(pos) = 1/A
  const float pdf_area = 1.0f/((maxa-mina)*(maxb-minb));
  return mf_mul(pdf, mf_set1(pdf_area));
}

mf_t shader_sky_sample_next_event(path_t *p)
{
  return rt.shader->skyshader.sample(p, rt.shader->skyshader.data);
}

mf_t shader_sky_pdf_next_event(const path_t *p, int v)
{
  return rt.shader->skyshader.pdf(p, v, rt.shader->skyshader.data);
}


mf_t shader_sky_eval(const path_t *p, int v)
{
  return rt.shader->skyshader.eval(p, v, rt.shader->skyshader.data);
}


int shader_shape_init(uint32_t shapeid, struct shader_so_t *self)
{
  if(self->shape_init) return self->shape_init(shapeid, self);
  return 0;
}

mf_t shader_pdf(const path_t *p, int v)
{
  int shader = p->v[v].hit.shader;
  if(shader < 0 && (p->v[v].mode & s_volume))
    shader = p->v[v].interior.shader;
  assert(shader >= 0);
  return rt.shader->shader[shader].pdf(p, v, v, v+1, rt.shader->shader[shader].data);
}

mf_t shader_pdf_adj(const path_t *p, int v)
{
  return rt.shader->shader[p->v[v].hit.shader].pdf(p, v+1, v, v, rt.shader->shader[p->v[v].hit.shader].data);
}

float shader_prepare(path_t *p, int v)
{
  if(p->v[v].flags & s_environment)
  {
    memset(&p->v[v].shading, 0, sizeof(vertex_shading_t));
    p->v[v].shading.roughness = 1.0f;
    p->v[v].shading.em = shader_sky_eval(p, v);
    return 1.0f;
  }
  // initialise tangent frames and diffgeo
  manifold_init(p, v);

  // no valid primid: in volumes or on sensor?
  if(primid_invalid(p->v[v].hit.prim))
  {
    // in participating medium,
    // no regular normal/primitive considerations hold here.
    
    // TODO: merge with below more elegantly.
    if(p->e[v].vol.shader >= 0)
    {
      p->v[v].interior = p->e[v].vol;
      if(rt.shader->shader[p->e[v].vol.shader].volume_enabled &&
         rt.shader->shader[p->e[v].vol.shader].volume_enabled(rt.shader->shader[p->e[v].vol.shader].data))
      { // heterogeneous
        memset(&p->v[v].shading, 0, sizeof(vertex_shading_t));
        p->v[v].shading.roughness = 1.0f;
        float res = 1.0f;
        shader_so_t *s = rt.shader->shader + p->e[v].vol.shader;
        if(s->prepare) res = fminf(res, s->prepare(p, v, s->data));
        // all other interactions (sample/pdf/eval) should be on the volume shader, too:
        p->v[v].hit.shader = p->e[v].vol.shader;
        return res;
      }
      else if(v > 0)
      { // homogeneous: copy over vertex volume data from previous vertex, if any
        // assume glossy lobe here (not pretty):
        p->v[v].material_modes = s_volume | s_glossy;
        return 1.0f;
      }
    }
    assert(p->v[v].mode & s_sensor);
    return 1.0f; // camera sensor
  }
  p->v[v].hit.shader = prims_shader(rt.prims, p->v[v].hit.prim);

  // default init:
  memset(&p->v[v].shading, 0, sizeof(vertex_shading_t));
  p->v[v].shading.roughness = 1.0f;

  // default init shape interior to vacuum
  path_volume_vacuum(&p->v[v].interior);

  // call shader's prepare func.
  float res = 1.0f;
  shader_so_t *s = rt.shader->shader + p->v[v].hit.shader;
  if(s->prepare) res = fminf(res, s->prepare(p, v, s->data));

  // pre-cache ior ratio for geometric manifolds:
  p->v[v].diffgeo.eta = path_eta_ratio(p, v);

#ifdef SHADER_ROUGHENING
  // horrible hack: use roughening after v[1] to regularise path space.
  // not consistent for bidirectional estimators.
  const float r = p->v[v].shading.roughness;
  const float t = CLAMP(v-1, 0, 100);
  const float mr = 1.7f;
  p->v[v].shading.roughness = CLAMP(r + (mr*t*t)/(100.*100.+1.), 0.f, 1.f);
#endif 
#if 0
  if(p->temperature > 0.f)
  { // melt measurement contribution function:
    const float r = p->v[v].shading.roughness;
    const float t = p->temperature;
    const float mr = 0.3f;
    p->v[v].shading.roughness = CLAMP(r + (mr*t*t)/(100.*100.+1.), 0.f, 1.f);
  }
#endif

  return res;
}

void shader_exterior_medium(path_t *p)
{
  if(rt.shader->exterior_medium_shader < 0)
  {
    path_volume_vacuum(&p->v[0].interior);
    p->e[0].vol = p->v[0].interior;
  }
  else
  {
    shader_so_t *s = rt.shader->shader + rt.shader->exterior_medium_shader;
    path_volume_vacuum(&p->v[0].interior);

    vertex_scattermode_t mode = p->v[0].mode;
    vertex_scattermode_t material_modes = p->v[0].material_modes;
    vertex_shading_t shading = p->v[0].shading;
    if(s->prepare) s->prepare(p, 0, s->data);
    p->v[0].mode = mode; // these would be overwritten by prepare()
    p->v[0].material_modes = material_modes;
    p->v[0].interior.shader = rt.shader->exterior_medium_shader;
    p->v[0].shading = shading;
    p->e[0].vol = p->v[0].interior;
  }
}

mf_t shader_brdf(path_t *p, int v)
{
  int shader = p->v[v].hit.shader;
  if(shader < 0 && (p->v[v].mode & s_volume))
    shader = p->v[v].interior.shader;
  assert(shader >= 0);
  return rt.shader->shader[shader].brdf(p, v, rt.shader->shader[shader].data);
}

mf_t shader_sample(path_t *p)
{
  const int v = p->length-1;
  assert(p->v[v].hit.shader >= 0);
  assert(p->v[v].hit.shader < rt.shader->num_shaders);
  const mf_t throughput = rt.shader->shader[p->v[v].hit.shader].sample(p, rt.shader->shader[p->v[v].hit.shader].data);
  normalise(p->e[v+1].omega); // super paranoid, but apparently we do get drift from higher bounces otherwise (cone intersection w/ lt shows it).
  // sidedness consistency check. we cannot let them pass or else our self intersection bias will fail, to say the least.
  const float dt = ((p->v[v].flags & s_inside) ? -1 : 1) * dotproduct(p->v[v].hit.gn, p->e[v+1].omega);
  if(((p->v[v].mode & s_reflect)  && (dt < 0.f)) ||
     ((p->v[v].mode & s_transmit) && (dt > 0.f)))
    return mf_set1(0.0f);
  return throughput;
}

int shader_inverse_sample(const path_t *p, int v, float *r_omega_x, float *r_omega_y, float *r_scatter_mode)
{
  assert(p->v[v].hit.shader >= 0);
  assert(p->v[v].hit.shader < rt.shader->num_shaders);
  if(rt.shader->shader[p->v[v].hit.shader].inverse_sample)
    return rt.shader->shader[p->v[v].hit.shader].inverse_sample(p, v, r_omega_x, r_omega_y, r_scatter_mode, rt.shader->shader[p->v[v].hit.shader].data);
  else
  { // not supported by this shader :(
    *r_omega_x = *r_omega_y = *r_scatter_mode = -1.0f;
    return 0;
  }
}

shader_t *shader_init(FILE* fd)
{
  shader_t *s = malloc(sizeof(*s));
  // :( but mult shader init needs it.
  rt.shader = s;

  s->num_handles = 0;
  s->exterior_medium_shader = -1;
  s->skyshader.eval = &sky_cloudy;
  s->skyshader.sample = sky_cloudy_sample;
  s->skyshader.pdf = &sky_cloudy_pdf;
  s->skyshader.black = 0;
  s->skyshader.sun = 0;
  s->skyshader.data = 0;
  s->skyshader.cleanup = 0;
  char name[1024];
  char obj[1024];
  int dreggn = 0;
  dreggn = fscanf(fd, "%s", name);
  if(dreggn == -1) fprintf(stderr, "gcc stinks\n");
  if(!strncmp(name, "daylight", 9))
  {
    s->skyshader.eval = &sky_daylight;
    s->skyshader.sample = &sky_daylight_sample;
    s->skyshader.pdf = &sky_daylight_pdf;
    s->skyshader.init = &sky_daylight_init;
    s->skyshader.init(fd, &(s->skyshader.data));
    s->skyshader.sun = 1;
  }
  else if(!strncmp(name, "black", 5))
  {
    s->skyshader.eval = &sky_black;
    s->skyshader.black = 1;
    dreggn = fscanf(fd, "%*[^\n]\n");
  }
  else if((!strncmp(name, "cloudy_sky", 10)) ||
          (!strncmp(name, "cloudy", 6)) ||
          (!strncmp(name, "clear_sky", 9)))
    dreggn = fscanf(fd, "%*[^\n]\n");
  else
  {
    sprintf(obj, "lib%s.so", name);
    s->handle[0] = dlopen(obj, RTLD_LAZY | RTLD_LOCAL);
    if(!dlerror())
    {
      s->skyshader.init = (init_t) dlsym(s->handle[0], "init");
      if(!dlerror()) s->skyshader.init(fd, &(s->skyshader.data));
      else dreggn = fscanf(fd, "%[^\n]\n", obj);
      s->skyshader.sample = (sky_sample_t) dlsym(s->handle[0], "sample");
      if(dlerror()) s->skyshader.sample = &shader_sky_sample_d;
      s->skyshader.pdf = (sky_pdf_t) dlsym(s->handle[0], "pdf");
      if(dlerror()) s->skyshader.pdf = &shader_sky_pdf_d;
      s->skyshader.cleanup = (cleanup_t) dlsym(s->handle[0], "cleanup");
      if(dlerror()) s->skyshader.cleanup = 0;
      s->skyshader.eval = (sky_eval_t) dlsym(s->handle[0], "eval");
      s->num_handles = 1;
      strncpy(s->dlname[0], name, 19);
      if(dlerror())
      {
        s->skyshader.eval = &sky_cloudy;
        s->skyshader.sample = &sky_cloudy_sample;
        s->skyshader.pdf = &sky_cloudy_pdf;
        dlclose(s->handle[0]);
        s->num_handles = 0;
      }
    }
    else
    {
      printf("[shader_init] failed to load sky shader `%s'\n", obj);
      dreggn = fscanf(fd, "%*[^\n]\n");
    }
  }
  s->num_shaders = 0;
  dreggn = fscanf(fd, "%d\n", &(s->num_shaders));
  s->shader = (shader_so_t *)malloc(sizeof(shader_so_t)*s->num_shaders);
  memset(s->shader, 0, sizeof(shader_so_t)*s->num_shaders);
  for(int k=0;k<s->num_shaders;k++)
  {
    s->shader[k].data = NULL;
    dreggn = fscanf(fd, "%s", name);
    int prepare_shader = 0;
    void* h = NULL;
    for(int i=0;i<s->num_handles;i++)
      if(strncmp(s->dlname[i], name, 19) == 0) h = s->handle[i];
    if(!h)
    {
      if(s->num_handles < 64)
      {
        sprintf(obj, "lib%s.so", name);
        h = s->handle[s->num_handles] = dlopen(obj, RTLD_LAZY | RTLD_LOCAL);
        // fprintf(stderr, "[shader_init] loading shader %s\n", name);
      }
      else fprintf(stderr, "[shader_init] limited to 64 different shader.so!\n");
      if(strcmp(name, "diffuse") == 0) goto diffuse;
      if(strcmp(name, "exterior") == 0)
      {
        dreggn  = fscanf(fd, "%d", &s->exterior_medium_shader);
        int light = 0;
        dreggn += fscanf(fd, "%d", &light);
        // dreggn += fscanf(f, "%*[^\n]\n");
        fprintf(stderr, "[shader_init] setting exterior medium to shader %d light sampling %d\n", s->exterior_medium_shader, light);
        if(light && s->exterior_medium_shader >= 0)
        {
          if(s->exterior_medium_shader >= k)
          {
            fprintf(stderr, "[shader_init] ERROR: exterior medium referring to a future medium shader!\n");
            goto error;
          }
          lights_init_volume_light(s->shader+s->exterior_medium_shader);
        }
        goto diffuse; // need to zero out all callbacks for cleanup later on
      }
      if (!h)
      {
        fprintf(stderr, "[shader_init] could not open %s in line %d.\n", obj, k+2);
        goto error;
      }
      dlerror();
      strncpy(s->dlname[s->num_handles], name, 19);
      s->num_handles++;
    }
    //printf("[shader_init] loading `%s'\n", name);
    s->shader[k].brdf = (brdf_t) dlsym(h, "brdf"); // surface shader, if brdf is given.
    if(dlerror()) { prepare_shader = 1; s->shader[k].brdf = &brdf_d; }
    if(!prepare_shader)
    {
      s->shader[k].sample = (sample_t) dlsym(h, "sample");
      if(dlerror()) { fprintf(stderr, "[shader_init] missing sample function in shader %s!\n", name); goto error; }
      s->shader[k].pdf = (pdf_t) dlsym(h, "pdf");
      if(dlerror()) { fprintf(stderr, "[shader_init] missing pdf function in shader %s!\n", name); goto error; }
    }
    s->shader[k].inverse_sample = (inverse_sample_t) dlsym(h, "inverse_sample");
    s->shader[k].prepare = (prepare_t) dlsym(h, "prepare");
    s->shader[k].shape_init = (shape_init_t) dlsym(h, "shape_init");

    s->shader[k].cleanup = (cleanup_t) dlsym(h, "cleanup");

    s->shader[k].volume_transmittance = dlsym(h, "volume_transmittance");
    s->shader[k].volume_sample = dlsym(h, "volume_sample");
    s->shader[k].volume_pdf = dlsym(h, "volume_pdf");
    s->shader[k].volume_pdf_adj = dlsym(h, "volume_pdf_adj");
    s->shader[k].volume_enabled = dlsym(h, "volume_enabled");

    s->shader[k].volume_sample_nee = dlsym(h, "volume_sample_nee");
    s->shader[k].volume_pdf_nee = dlsym(h, "volume_pdf_nee");
    s->shader[k].volume_sample_light = dlsym(h, "volume_sample_light");
    s->shader[k].volume_pdf_light = dlsym(h, "volume_pdf_light");
    s->shader[k].volume_sample_fnee_direction = dlsym(h, "volume_sample_fnee_direction");
    s->shader[k].volume_pdf_fnee_direction = dlsym(h, "volume_pdf_fnee_direction");
    s->shader[k].volume_pdf_fnee_direction_adj = dlsym(h, "volume_pdf_direction_adj");

    s->shader[k].init = (init_t) dlsym(h, "init");
    if(dlerror()) s->shader[k].init = &init_d;
    if(s->shader[k].init(fd, &(s->shader[k].data))) goto error;

    continue;
error:
    printf("[shader_init] WARN: can't open shader: `%s' in line %d\n", name, k+3);
    printf("  %s\n", dlerror());
diffuse:
    // clear input line:
    init_d(fd, &(s->shader[k].data));
    s->shader[k].init = 0;
    s->shader[k].shape_init = 0;
    s->shader[k].prepare = &prepare_d;
    s->shader[k].cleanup = 0;
    s->shader[k].pdf = &pdf_d;
    s->shader[k].sample = &sample_d;
    s->shader[k].inverse_sample = 0;
    s->shader[k].brdf = &brdf_d;
    s->shader[k].data = s;

    s->shader[k].volume_pdf = 0;
    s->shader[k].volume_pdf_adj = 0;
    s->shader[k].volume_sample = 0;
    s->shader[k].volume_transmittance = 0;
    s->shader[k].volume_sample_nee = 0;
    s->shader[k].volume_pdf_nee = 0;
    s->shader[k].volume_sample_light = 0;
    s->shader[k].volume_pdf_light = 0;
    s->shader[k].volume_sample_fnee_direction = 0;
    s->shader[k].volume_pdf_fnee_direction = 0;
    s->shader[k].volume_pdf_fnee_direction_adj = 0;
  }

#if 0 // sorry can't do this (volume_enabled may access negative indices for mult shaders)
  // consistency check on heterogenous volumes:
  for(int k=0;k<s->num_shaders;k++)
  {
    if(s->shader[k].volume_enabled && s->shader[k].volume_enabled(s->shader[k].data))
    {
      const int vol_check =
        ((s->shader[k].volume_transmittance != 0)&1) +
        ((s->shader[k].volume_sample != 0)&1) +
        ((s->shader[k].volume_pdf != 0)&1) +
        ((s->shader[k].volume_pdf_adj != 0)&1);

      if(vol_check < 4)
      {
        fprintf(stderr, "[shader_init] WARN: shader %d does not define all of the 4 volume scattering functions!\n", k);
        fprintf(stderr, "[shader_init] WARN: disabling heterogeneous scattering!\n");
        goto no_volume;
      }
    }
    else
    { // no volume enabled
no_volume:
      s->shader[k].volume_transmittance = 0;
      s->shader[k].volume_sample = 0;
      s->shader[k].volume_pdf = 0;
      s->shader[k].volume_pdf_adj = 0;
    }
  }
#endif
  return s;
}

void shader_cleanup(shader_t *s)
{
  for(int k=0;k<s->num_shaders;k++)
    if(s->shader[k].cleanup) s->shader[k].cleanup(s->shader[k].data);
    else if(s->shader[k].init && s->shader[k].data) free(s->shader[k].data);
  free(s->shader);
  if(s->skyshader.cleanup)
    s->skyshader.cleanup(s->skyshader.data);
  else if(s->skyshader.data) free(s->skyshader.data);
  for(int k=0;k<s->num_handles;k++) dlclose(s->handle[k]);
  free(s);
}

