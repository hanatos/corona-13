#pragma once

#include "pathspace.h"
#include "pointsampler.h"
#include "sampler_common.h"
#include "shader.h"
#include "prims.h"
#include <stdint.h>

typedef struct lights_t lights_t;

// global initialisation
lights_t *lights_init();

// ..and cleanup
void lights_cleanup(lights_t *s);

// attach a volume emitter to the global list
void lights_init_volume_light(shader_so_t *vol);

// initialise the give shapeid as light source, use L as average emitted
// radiance for importance sampling
void lights_init_light(const uint32_t shapeid, const float L);

// called once before every progression
void lights_prepare_frame();

// returns the probability to sample environment, geo, or volume lights
// as next event estimation sample from vertex p->v[v].
// if v is the end point of a path and mode & s_emit, it is interpreted
// as the probability to start a light ray instead of nee.
void lights_pdf_type(const path_t *p, const int v, mf_t p_egv[3]);

// computes the spectral pdf of next event estimation for p->v[v]
mf_t lights_pdf_next_event(const path_t *p, int v);

// returns throughput into direction of segment, no visibility is checked here.
// will construct vertex p->v[p->length] and not increment length
mf_t lights_sample_next_event(path_t *p);

// sample a light ray at p->v[0] and p->e[1]
mf_t lights_sample(path_t *p);

// returns the pdf of sampling the first or last segment
// of the path (depending on v) via lights_sample()
mf_t lights_pdf(const path_t *p, int v);

// evaluate emitted radiance of given path vertex v,
// in direction towards the sensor.
mf_t lights_eval_vertex(path_t *path, int v);

// volume sampling wrappers
mf_t light_volume_pdf_nee(const path_t *p, int v);
mf_t light_volume_pdf_fnee_direction(const path_t *p, int v);
mf_t light_volume_pdf(const path_t *p, int v);
mf_t light_volume_sample_nee(path_t *p, int v);
void light_volume_sample_fnee_direction(path_t *p, int v);



// evaluate emitted radiance of light source at the end point.
// auto-detect light and path tracing cases.
static inline mf_t lights_eval(path_t *path)
{
  int v = path->length-1;
  if(path->v[0].mode & s_emit) v = 0;
  return lights_eval_vertex(path, v);
}

// sample all types of light sources (envmap, volume, geo)
static inline mf_t light_sampler_sample(path_t* path)
{
  mf_t p_egv[3];
  lights_pdf_type(path, 0, p_egv);
  const float rand = pointsampler(path, s_dim_envmapvsarea);
  if(rand < mf(p_egv[0],0))
  {
    path->v[0].flags = s_environment;
    path->v[0].throughput = mf_div(shader_sky_sample(path), mf_set1(mf(p_egv[0], 0)));
    path->v[1].pdf = mf_mul(path->v[1].pdf, p_egv[0]);
  }
  else if(rand < mf(p_egv[0],0)+mf(p_egv[1],0))
  {
    path->v[0].throughput = mf_div(lights_sample(path), mf_set1(mf(p_egv[1],0)));
    path->v[1].pdf = mf_mul(path->v[1].pdf, p_egv[1]);
  }
  else if(rand < mf(p_egv[0],0)+mf(p_egv[1],0)+mf(p_egv[2],0))
  {
    assert(0 && "wire light tracing for volumes!");
    // path->v[0].throughput = rt.lights->vol->lights_sample(path)/rt.lights->p_vol;
    // path->v[1].pdf *= rt.lights->p_vol;
  }
  else
  {
    fprintf(stderr, "[lights] broken random number %g\n", rand);
  }
  return path->v[0].throughput;
}

// returns the pdf of sampling a light path segment at v[0]
// considering all three options (env geo vol)
static inline mf_t light_sampler_pdf_extend(const path_t* path)
{
  mf_t pdf = mf_set1(1.0f);
  mf_t p_egv[3];
  lights_pdf_type(path, 0, p_egv);
  if(path->v[0].flags & s_environment)
    pdf = mf_mul(p_egv[0], shader_sky_pdf(path, 1));
  else if(primid_invalid(path->v[0].hit.prim))
    // XXX fix this for light tracing!
    assert(0);
    // pdf = rt.lights->p_vol*rt.lights->vol->volume_pdf_nee(path, 0, rt.lights->vol->data);
  else // geo light
    pdf = mf_mul(p_egv[1], lights_pdf(path, 0));
  return pdf;
}

// returns the pdf of sampling a light path segment from behind (v+1..v)
static inline mf_t light_sampler_pdf_extend_adjoint(const path_t *path, int v)
{
  mf_t p_egv[3];
  lights_pdf_type(path, 0, p_egv);
  if(path->v[v+1].flags & s_environment)
    return mf_mul(p_egv[0], shader_sky_pdf(path, v+1));
  else if(primid_invalid(path->v[v+1].hit.prim) && (path->v[v+1].mode & s_emit))
    assert(0 && "wire light tracing stuff");
    // return = rt.lights->p_vol*rt.lights->vol->volume_pdf_nee(path, v, rt.lights->vol->data);
  else if(path->v[v+1].mode & s_emit)
    return mf_mul(p_egv[1], lights_pdf(path, v+1));
  return mf_set1(0.0f);
}

