#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"

#include <math.h>

#if 0
static inline float t_rand(uint64_t *state)
{
  uint64_t s1 = state[0];
  uint64_t s0 = state[1];
  state[0] = s0;
  s1 ^= s1 << 23;
  s1 ^= s1 >> 17;
  s1 ^= s0;
  s1 ^= s0 >> 26;
  state[1] = s1;
  // return (state0 + state1) / ((double)((uint64_t)-1) + 1.0);
  uint32_t v = 0x3f800000 | ((state[0]+state[1])>>41); // faster than double version.
  return (*(float*)&v) - 1.0f;
}
#endif

static inline int
mvnee_possible(const path_t *p, const int v)
{
  // if there are no non-specular components, fail next event estimation.
  if(!(p->v[v].material_modes & (s_diffuse | s_glossy)))
    return 0;
  if(v == 0) return 0; // as nee: do not connect light source directly (though we could. should we?)
  if(v + 2 >= PATHSPACE_MAX_VERTS) return 0; // not enough vertices on list
  // TODO: can we test whether the vertex will potentially be in a volume already?
  // v[v].interior will potentially be inside the other object if we transmit,
  // but we could probably bail out if neither e[v].vol nor v[v].interior have a volume.
  return 1;
}

// this returns the product vertex area measure pdf of sampling
// both end vertices with mvnee.
static inline mf_t
mvnee_pdf(const path_t* p, int v)
{
  if(p->length < 3) return mf_set1(0.0f); // no next event for 2-vertex paths.
  assert(v == 0 || v == p->length-1);
  int v0 = v ? v-2 : v+2;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int v1 = v ? v-1 : v+1;  // middle vertex created by mvnee
  int e0 = v ? v-1 : 2;    // towards double scattering vertex
  int e1 = v ? v   : 1;    // towards last point
  if(!mvnee_possible(p, v0)) return mf_set1(0.0f);
  if(!primid_invalid(p->v[v1].hit.prim)) return mf_set1(0.0f);

  mf_t p_egv[3];
  lights_pdf_type(p, v0, p_egv);

  mf_t res = mf_set1(0.0f);
  float cos_theta = dotproduct(p->e[e0].omega, p->e[e1].omega);
  if(cos_theta < 0.0f) return mf_set1(0.0f);

  // light tracer connects to camera
  if(((p->v[0].mode & s_emit)>0) ^ (v==0))
    res = view_cam_pdf_connect(p, v);
  else if(p->v[v].flags & s_environment)
    res = mf_mul(p_egv[0], shader_sky_pdf_next_event(p, v));
  else if(primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf_any(mf_gt(p->v[v].shading.em, mf_set1(0.0f))) && mf(p_egv[2], 0) > 0)
    res = mf_mul(p_egv[2], light_volume_pdf_nee(p, v));
  else if(!primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf(p_egv[1], 0) > 0)
    res = mf_mul(p_egv[1], lights_pdf_next_event(p, v));
  else return mf_set1(0.0f);

  // res is now the vertex area pdf of sampling the light source position v[v]
  const float distv[] = {
    p->v[v].hit.x[0] - p->v[v0].hit.x[0],
    p->v[v].hit.x[1] - p->v[v0].hit.x[1],
    p->v[v].hit.x[2] - p->v[v0].hit.x[2]};
  const double s = sqrt(dotproduct(distv, distv));
  const double sinh = MAX(1e-8, 1.0 - cos_theta*cos_theta);
  const double theta = acos(CLAMP(cos_theta, 0.0, 1.0));
  const float g = p->v[v1].interior.mean_cos;
  const double hg_pdf = sample_eval_hg_fwd(g, p->e[e0].omega, p->e[e1].omega);
  const double sinc = theta < 1e-7 ? 1.0 : sinh / theta;
  return mf_mul(res, // v[v], the nee vertex
      mf_set1(fabs(hg_pdf/(p->e[e0].dist*p->e[e0].dist*p->e[e1].dist*p->e[e1].dist)
          * s * sinc)));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// next event estimation of the reverse path from v+1
static inline mf_t
mvnee_pdf_adjoint(const path_t *path, int v)
{
  // luckily these guys know already about v==0 or v==length-1
  return mvnee_pdf(path, v);
}

// sample next event. returns != 0 on failure
static inline int
mvnee_sample(path_t *p)
{
  assert(p->length >= 2); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  if(p->length >= PATHSPACE_MAX_VERTS-1) return 1; // need to append two new vertices
  const int v = p->length; // constructing new vertex here (and at v+1)

  if(!mvnee_possible(p, v-1)) goto fail;

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
  p->v[v].tech = s_tech_mvnee;

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
    edf = light_volume_sample_nee(p, v);
    if(!mf_any(mf_gt(edf, mf_set1(.0f)))) goto fail;
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[2]);
    edf = mf_div(edf, mf_set1(mf(p_egv[2],0)));
    // init segment
    for(int k=0;k<3;k++)
      p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
    p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
    for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
  }

  // ask edf and bsdf for their consent:
  if(!mf_any(mf_gt(edf, mf_set1(0.0f)))) goto fail;
  // compute brdf and throughput:
  mf_t bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1
  if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.

  // determine side of surface and volume from that (brdf sets mode)
  if(path_edge_init_volume(p, v)) goto fail;

  // move vertex v to v+1, free up space for v[v], our in-between vertex.
  p->e[v+1] = p->e[v];
  p->v[v+1] = p->v[v];
  // get volume properties:
  if(p->e[v].vol.shader == -1) goto fail; // no volume no mvnee.
  const float g = p->e[v].vol.mean_cos;
  // fprintf(stderr, "mean cos %g\n", g);
  // sample new vertex v on edge between the two:
  p->length = v+1; // instruct pointsampler to get new dimensions
  float hg_r1 = pointsampler(p, s_dim_nee_x); // used for cosh
  float hg_r2 = pointsampler(p, s_dim_nee_y); // usef for isotropic phi
  // sample hg with g
  float hg_pdf = 0.0f, omega[3];
  sample_hg_fwd(g, hg_r1, hg_r2, omega, &hg_pdf);
  const double cosh = omega[0];
  const double sinh = sqrt(MAX(0.0, 1.0-cosh*cosh));
  const double theta = acos(CLAMP(omega[0], 0.0, 1.0));

  // the cosine from phase function sampling defines the angle between e[v].w and e[v+1].w.
  const float dist_r = pointsampler(p, s_dim_nee_light1);
  // sample our new inverse cdf:
  const float t = CLAMP(cosf(theta - dist_r * theta) * sinf(dist_r * theta) / MAX(1e-5, sinh), 0.0f, 1.0f);
  // const float t = dist_r; // uniformly

  // angle around axis of symmetry
  const float phi = 2.0f*M_PI*hg_r2;
  const float sin_phi= sinf(phi);
  const float cos_phi= cosf(phi);
  const float sinh2 = 1.0-omega[0]*omega[0];
  const double r = sqrt(1.0/(4.0*sinh2) - (0.5-t)*(0.5-t)) - sqrt(1.0/(4.0*sinh2) - 0.25);
  const double s = p->e[v].dist;

  // initialise edges and adjust p->v[v].hit.x accordingly.
  // create onb around e->omega (straight connection)
  float onb_u[3], onb_v[3];
  get_onb(p->e[v].omega, onb_u, onb_v);
  for(int k=0;k<3;k++)
    p->v[v].hit.x[k] = p->v[v-1].hit.x[k] * (1.0f-t) + p->v[v+1].hit.x[k] * t + 
      s*r*(sin_phi * onb_u[k] + cos_phi * onb_v[k]);
  // also make it a point in the volume, shader_prepare will pick up the rest below
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = p->e[v].vol.shader;

  // recreate omega/dist on the edges
  for(int k=0;k<3;k++)
    p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
  p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
  for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
  for(int k=0;k<3;k++)
    p->e[v+1].omega[k] = p->v[v+1].hit.x[k] - p->v[v].hit.x[k];
  p->e[v+1].dist = sqrtf(dotproduct(p->e[v+1].omega, p->e[v+1].omega));
  for(int k=0;k<3;k++) p->e[v+1].omega[k] *= 1./p->e[v+1].dist;

  // test visibility of both new segments:
  if(!path_visible(p, v  )) goto fail;
  if(!path_visible(p, v+1)) goto fail;

  bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1 according to new direction. also sets mode.
  if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.
  if(path_edge_init_volume(p, v))   goto fail; // init volumes according to new mode
  shader_prepare(p, v);
  mf_t phasef = shader_brdf(p, v); // set mode
  (void)phasef;
  if(path_edge_init_volume(p, v+1)) goto fail;

  shader_prepare(p, v+1);
  // sampled point so far at the rims of the volume that we fell
  // off it trying to determine the emission
  if((p->v[v].mode & s_volume) && mf_all(mf_eq(p->v[v].interior.mu_t, mf_set1(0.0f))))
    goto fail;

  // XXX hack works for lt only: reconnect camera:
  p->length = v+2; // let cam know what's the last vertex
  // we sampled this already during nee. in particular, don't use any more
  // random numbers on the last vertex (there aren't any, nee has to be
  // performed one vertex back)
  p->sensor.aperture_set = 1;
  if(p->v[0].mode & s_emit)
    edf = view_cam_connect(p);
  // TODO: also reconnect to light (in case of cosine power EDF)

  // compute transmittance and egde contribution
  const mf_t transmittance =
    shader_vol_transmittance(p, v+1) *
    shader_vol_transmittance(p, v);

  // TODO: multiply by shader_brdf[v] / pdf_hg to account for backward scattering normalisation

  // these are the terms we don't sample:
  const md_t f = md_mul(md_mul(md_mul(mf_2d(bsdf), mf_2d(transmittance)), md_mul(mf_2d(edf), mf_2d(p->v[v].interior.mu_s))),
    md_set1(path_lambert(p, v-1, p->e[v].omega)*path_lambert(p, v+1, p->e[v+1].omega)));
  // this is when just sampling the phase function, the remaining jacobian:
  // const md_t throughput = md_mul(f, md_set1(GJ));
  // this is the remaining weight when sampling both theta and t:
  const double isinc = theta < 1e-7 ? 1.0 : theta / sinh;
  const md_t throughput = md_mul(f, md_set1(isinc/s));
  p->v[v].throughput = md_2f(md_mul(throughput, mf_2d(p->v[v-1].throughput)));

  p->v[v+1].throughput = p->v[v].throughput; // just set both to the same thing
  p->throughput = p->v[v+1].throughput; // already contains light/sensor edf


#if 0 // fireflies seem to come from paths like:
 LG [0] throughput 39982840.000000 eta 1.000000
    dist 598.446411 ior 1.000000
 TG [1] throughput 39982840.000000 eta 0.750873
 |  tau 0.786251 em 0.000000 pdf 0.009621 dist 19.652620 ior 1.331783
 oG [2] throughput 3266103296.000000 eta 1.000000
 |  tau 0.544346 em 0.000000 pdf 1.000000 dist 49.701405 ior 1.331783
 oG [3] throughput 9.857325 eta 1.000000
 |  tau 0.096501 em 0.000000 pdf 1.000000 dist 191.084076 ior 1.331783
 E  [4] throughput 9.857325 eta 1.00000
 where oG [2] happens to have an angular configuration that leaves too much of the phase function :(
  if(throughput > 1e-9)
  { // XXX DEBUG
    path_print(p, stdout);
    goto fail;
  }
#endif

  if(0)
  {
fail:
    p->v[v].pdf = mf_set1(0.0f);
    p->throughput = mf_set1(0.0f);
    p->v[v].throughput = mf_set1(0.0f);
    p->v[v].flags = s_none;
    p->v[v].mode = s_absorb;
    p->v[v].rand_cnt = s_dim_num_nee;
    p->v[v+1] = p->v[v];
    // need two more, because pop will also pop two
    p->length = v+2; // constructed vertex v (even if it absorbs)
    return 0;
  }
  p->v[v  ].rand_cnt = s_dim_num_nee;
  p->v[v+1].rand_cnt = s_dim_num_nee;
  p->length = v+2; // constructed vertex v and v+1

  // vertex area pdf of sampling the light source position:
  const mf_t pdf_nee = p->v[v+1].pdf;
  p->v[v+1].total_throughput = p->v[v+1].throughput; // store for path_pop()
  p->v[v+1].pdf = pdf_nee;
  
  // compute vertex area measure pdf of double scatter vertex v[v]:
  const double sinc = theta < 1e-7 ? 1.0 : sinh / theta;
  p->v[v].pdf = mf_set1(hg_pdf * fabs(1.0/(p->e[v].dist*p->e[v].dist*p->e[v+1].dist*p->e[v+1].dist) *
        s * sinc));

  return 0;
}
