#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"

static inline int nee_possible(const path_t *p, const int v)
{
  // if there are no non-specular components, fail next event estimation.
  if(!(p->v[v].material_modes & (s_diffuse | s_glossy)))
    return 0;
  // don't connect from inside closed shapes without enclosed light sources!
  if((p->v[v].hit.shader != rt.shader->exterior_medium_shader) &&
      primid_invalid(p->v[v].hit.prim)) return 0;
  return 1;
}

static inline mf_t nee_pdf_nee(const path_t* p, int v)
{
  if(p->length < 3) return mf_set1(0.0f); // no next event for 2-vertex paths.
  assert(v == 0 || v == p->length-1);
  int v1 = v ? v-1 : v+1;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  if(!nee_possible(p, v1)) return mf_set1(0.0f);

  mf_t p_egv[3];
  lights_pdf_type(p, v1, p_egv);

  // light tracer connects to camera
  if(((p->v[0].mode & s_emit)>0) ^ (v==0))
    return view_cam_pdf_connect(p, v);
  else if(p->v[v].flags & s_environment)
    return mf_mul(p_egv[0], shader_sky_pdf_next_event(p, v));
  else if(primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf_any(mf_gt(p->v[v].shading.em, mf_set1(0.0f))) && mf(p_egv[2], 0) > 0)
#ifndef FNEE
    return mf_mul(p_egv[2], light_volume_pdf_nee(p, v));
#else
    return mf_set1(0.0f);
#endif
  // XXX TODO: check whether we're in trouble because s_emit is also set on vertices which are /not/ emissive but the incoming segment is!
  else if(!primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf(p_egv[1], 0) > 0)
    return mf_mul(p_egv[1], lights_pdf_next_event(p, v));
  return mf_set1(0.0f);
}

static inline mf_t nee_pdf_fnee(const path_t* p, int v)
{
#ifndef FNEE
  return mf_set1(0.0f);
#else
  mf_t p_egv[3];
  lights_pdf_type(p, v1, p_egv);
  if(p->length < 3 || !(mf(p_egv[2], 0) > 0.f)) return mf_set1(0.0f); // no next event for 2-vertex paths.
  assert(v == 0 || v == p->length-1);
  int v1 = v ? v-1 : v+1;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int e  = v ? v : v+1;    // edge leading up to the next event vertex, again v==0 is adjoint
  if(!nee_possible(p, v1)) return mf_set1(0.0f);

  if (p->e[e].vol.shader == -1)
    return mf_set1(0.0f);

  const float cos_theta = path_lambert(p, v-1, p->e[v].omega);
  mf_t pdf_dir = light_volume_pdf_fnee_direction(p, v);

  // this should be correct:
  const mf_t pdf_e = light_volume_pdf(p, v);
  return mf_mul(mf_mul(p_egv[2], mf_set1(1.f/cos_theta * path_G(p, v))), mf_mul(pdf_dir, pdf_e));
#endif
}

static inline mf_t nee_pdf(const path_t *path, int v)
{
  return mf_add(nee_pdf_nee(path, v), nee_pdf_fnee(path, v));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// next event estimation of the reverse path from v+1
static inline mf_t nee_pdf_adjoint(const path_t *path, int v)
{
  // luckily these guys know already about v==0 or v==length-1
  return nee_pdf(path, v);
}

// sample next event. returns != 0 on failure
static inline int nee_sample(path_t *p)
{
  assert(p->length >= 2); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  if(p->length >= PATHSPACE_MAX_VERTS) return 1;
  const int v = p->length; // constructing new vertex here

  if(!nee_possible(p, v-1)) goto fail;

  // sample vertex v as endpoint (lights or camera)
  mf_t edf = mf_set1(0.0f);
  mf_t p_egv[3];
  lights_pdf_type(p, v-1, p_egv);

  // constructing v[v] here:
  memset(p->v+v, 0, sizeof(vertex_t));
  memset(p->e+v, 0, sizeof(edge_t));
  p->v[v].pdf_mnee = mf_set1(-1.f);

  // instruct kelemen mlt to use new random numbers:
  p->v[v].rand_beg = p->v[v-1].rand_beg + p->v[v-1].rand_cnt;

  // set technique to next event estimation
  p->v[v].tech = s_tech_nee;

  const float rand = pointsampler(p, s_dim_nee_light1);
  if(p->v[0].mode & s_emit)
  { // connect to the camera
    edf = view_cam_connect(p);
  }
  else if(rand < mf(p_egv[0],0))
  { // connect to the envmap
    edf = shader_sky_sample_next_event(p);
    edf = mf_div(edf, mf_set1(mf(p_egv[0], 0)));
    // pdf is in solid angle
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[0]);
  }
  else if(rand < mf(p_egv[0],0)+mf(p_egv[1],0))
  { // connect to geo lights
    edf = lights_sample_next_event(p);
    // compensate for envmap sampling probability
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[1]);
    edf = mf_div(edf, mf_set1(mf(p_egv[1], 0)));
  }
  else if(mf(p_egv[2],0) > 0)
  { // connect to volume lights
#ifndef FNEE
    edf = light_volume_sample_nee(p, v);
    if(!mf_any(mf_gt(edf, mf_set1(.0f)))) goto fail;
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[2]);
    edf = mf_div(edf, mf_set1(mf(p_egv[2],0)));
    // init segment
    for(int k=0;k<3;k++)
      p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
    p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
    for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
#else // volume forward sampling:
    light_volume_sample_fnee_direction(p, v);
    if(!mf_any(mf_gt(p->v[v].pdf, mf_set1(0.0f)))) goto fail;
    
    const mf_t bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1, and any extra ray epsilons if required.
    if(mf_all(mf_lte(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.

    const float cos_theta = path_lambert(p, v-1, p->e[v].omega);
    p->v[v].pdf = mf_mul(p->v[v].pdf, mf_mul(mf_set1(1.f/cos_theta), p_egv[2]));
    p->v[v].throughput = mf_mul(p->v[v-1].throughput, mf_div(bsdf,
          mf_mul(p->v[v].pdf, mf_set1(1.f/cos_theta * mf(p_egv[2],0)))));
    
    if(path_propagate(p, v, s_propagate_sample)) goto fail;

    // sampled point so far at the rims of the volume that we fell
    // off it trying to determine the emission
    if((p->v[v].mode & s_volume) && mf_all(mf_eq(p->v[v].interior.mu_t, mf_set1(0.0f))))
      goto fail;

    p->v[v].rand_cnt = s_dim_num_nee;
    p->length++; // constructed vertex v (even if it absorbs)
    
    p->v[v].pdf = mf_mul(p->v[v].pdf, mf_set1(path_G(p, v)));
    path_update_throughput(p, v);
    
    // now we can compute internal fnee-nee MIS weights
    // we do that by comparing projected solid angle (on geo lights)
    // or solid angle (on the envmap) pdfs. the fnee estimator gives
    // the correct answer when considering all paths, some of which ending
    // inside the volume. the nee estimator gives the fully correct result,
    // too, but never ends inside the volume. we thus need to weight the nee
    // estimate against all fnee paths, integrating out the transmittance/edge pdf:
    const mf_t pdf_fnee = p->v[v].pdf;
    const mf_t pdf_nee = nee_pdf_nee(p, v);
    const mf_t weight = mf_div(pdf_fnee, mf_fma(pdf_nee, p->e[v].pdf, pdf_fnee));

    p->throughput = mf_mul(p->throughput, weight);
    p->v[v].throughput = mf_mul(p->v[v].throughput, weight);

    p->v[v].total_throughput = p->v[v].throughput; // store for path_pop()
    p->v[v].pdf = mf_add(pdf_fnee, pdf_nee);
    return 0;
#endif
  }

  if(!mf_any(mf_gt(edf, mf_set1(0.0f)))) goto fail;
  // need to check visibility to new vertex, compute brdf and throughput:
  const mf_t bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1
  if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.

  // determine side of surface and volume from that (brdf sets mode)
  if(path_edge_init_volume(p, v)) goto fail;

  if(!path_visible(p, v)) goto fail;

  shader_prepare(p, v);
  // sampled point so far at the rims of the volume that we fell
  // off it trying to determine the emission
  if((p->v[v].mode & s_volume) && mf_all(mf_eq(p->v[v].interior.mu_t, mf_set1(0.0f))))
    goto fail;

  const float G = path_G(p, v);

  // compute transmittance and egde contribution
  const mf_t transmittance = shader_vol_transmittance(p, v);

  p->v[v].throughput = mf_mul(mf_mul(mf_mul(p->v[v-1].throughput, bsdf), mf_mul(transmittance, edf)), mf_set1(G));
  p->v[v].throughput = mf_add(p->v[v].throughput, mf_mul(mf_mul(p->v[v-1].throughput, bsdf), mf_div(mf_mul(p->e[v].contribution, mf_set1(G)), p->v[v].pdf)));

  p->throughput = p->v[v].throughput; // already contains light/sensor edf

  if(0)
  {
fail:
    p->v[v].pdf = mf_set1(0.0f);
    p->throughput = mf_set1(0.0f);
    p->v[v].throughput = mf_set1(0.0f);
    p->v[v].flags = s_none;
    p->v[v].mode = s_absorb;
    p->v[v].rand_cnt = s_dim_num_nee;
    p->length++; // constructed vertex v (even if it absorbs)
    return 0;
  }
  p->v[v].rand_cnt = s_dim_num_nee;
  p->length++; // constructed vertex v

  const mf_t pdf_nee = p->v[v].pdf;
  const mf_t pdf_fnee = nee_pdf_fnee(p, v);
  // interestingly these are both vertex area measure on geo light sources or
  // solid angle measure on the envmap:
  const mf_t weight = mf_div(pdf_nee, mf_add(pdf_nee, mf_div(pdf_fnee, transmittance)));
  p->throughput = mf_mul(p->throughput, weight);
  p->v[v].throughput = mf_mul(p->v[v].throughput, weight);
  p->v[v].total_throughput = p->v[v].throughput; // store for path_pop()
  // const float pdf2 = nee_pdf(p, v);
  // const float total_pdf = pdf_nee + pdf_fnee;
  // if(fabsf(pdf2 - total_pdf) > 1e-3f*MAX(pdf2, total_pdf))
  //   fprintf(stderr, "pdf lies: %g %g\n", pdf2, total_pdf);
  p->v[v].pdf = mf_add(pdf_nee, pdf_fnee);
  return 0;
}
