#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"

#include <math.h>


static inline int
equiangular_possible(const path_t *p, const int v)
{
  // if there are no non-specular components, fail next event estimation.
  if(!(p->v[v].material_modes & (s_diffuse | s_glossy)))
    return 0;
  // TODO: can we test whether the vertex will potentially be in a volume already?
  // v[v].interior will potentially be inside the other object if we transmit,
  // but we could probably bail out if neither e[v].vol nor v[v].interior have a volume.
  return 1;
}

// this returns the product vertex area measure pdf of sampling
// both end vertices with mvnee.
static inline mf_t
equiangular_pdf(const path_t* p, int v)
{
  if(p->length < 3) return mf_set1(0.0f); // no next event for 2-vertex paths.
  assert(v == 0 || v == p->length-1);
  int v0 = v ? v-2 : v+2;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int v1 = v ? v-1 : v+1;  // middle vertex created by mvnee
  int e0 = v ? v-1 : 2;    // towards double scattering vertex
  int e1 = v ? v   : 1;    // towards last point
  if(!equiangular_possible(p, v0)) return mf_set1(0.0f);
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

  double equi_pdf;
  float tolight[3];
  for(int k=0;k<3;k++) tolight[k] = p->v[v].hit.x[k] - p->v[v0].hit.x[k];
  const double delta = dotproduct(p->e[v].omega, tolight);
  float ortho[3];
  for(int k=0;k<3;k++)
    ortho[k] = p->v[v0].hit.x[k] + p->e[e0].omega[k] * delta - p->v[v].hit.x[k];
  const double D = sqrtf(dotproduct(ortho, ortho));

  double a = 0.0 - delta;
  double b = 10000.0 - delta; // sorry, this is an infinite volume :(

  const double t = p->e[e0].dist - delta;
  if(D > 1e-8)
  {
    double thetaA = atan(a/D);
    double thetaB = atan(b/D);
    equi_pdf = D/fabs(thetaA - thetaB)/(D*D + t*t);
  }
  else
  {
    equi_pdf = a*b / (b - a) / (t*t);
  }
  mf_t p_bsdf = shader_pdf(p, v1);

  return mf_mul(mf_mul(p_bsdf, res), mf_set1(equi_pdf * path_G(p, e0)));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// next event estimation of the reverse path from v+1
static inline mf_t
equiangular_pdf_adjoint(const path_t *path, int v)
{
  // luckily these guys know already about v==0 or v==length-1
  return equiangular_pdf(path, v);
}

// sample next event. returns != 0 on failure
static inline int
equiangular_sample(path_t *p)
{
  assert(p->length >= 2); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  if(p->length >= PATHSPACE_MAX_VERTS-1) return 1; // need to append two new vertices
  const int v = p->length; // constructing new vertex here (and at v+1)

  if(!equiangular_possible(p, v-1)) goto fail;

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

  // set technique
  p->v[v].tech = s_tech_equiangular;

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

  // move vertex v to v+1, free up space for v[v], our in-between vertex.
  p->e[v+1] = p->e[v];
  p->v[v+1] = p->v[v];
  // sample new vertex v on edge between the two:

  // ================================================================================================
  // sample bsdf, proceed as path_extend() but don't use distance sampling by transmittance:
  p->length = v;
  memset(p->v+v, 0, sizeof(vertex_t));
  memset(p->e+v, 0, sizeof(edge_t));

  // remember our random number offset
  p->v[v].rand_beg = p->v[v-1].rand_beg + p->v[v-1].rand_cnt;
  p->v[v].pdf = mf_set1(1.0f);
  p->v[v].throughput = p->v[v-1].throughput;

  // sample the bsdf at vertex v-1. also inits e[v] with volume information.
  p->v[v].throughput = mf_mul(p->v[v].throughput, shader_sample(p));
  if(path_edge_init_volume(p, v)) goto fail;
  if(p->e[v].vol.shader == -1) goto fail; // no volume no equiangular

  // this includes s_dim_russian_r, in case it is needed later
  // by means of calling path_russian_roulette(.).
  p->v[v].rand_cnt = s_dim_num_extend;

  // bsdf sampling failed?
  if(mf_all(mf_lte(p->v[v].throughput, mf_set1(0.0f)))) goto fail;

  // set distance by equiangular sampling!
  // stolen from implementation in smallpt: https://github.com/sakanaman/equi-angular-sampling.git
  float tolight[3];
  for(int k=0;k<3;k++) tolight[k] = p->v[v+1].hit.x[k] - p->v[v-1].hit.x[k];
  const double delta = dotproduct(p->e[v].omega, tolight);
  float ortho[3];
  for(int k=0;k<3;k++)
    ortho[k] = p->v[v-1].hit.x[k] + p->e[v].omega[k] * delta - p->v[v+1].hit.x[k];
  const double D = sqrtf(dotproduct(ortho, ortho));

  double a = 0.0 - delta;
  double b = 10000.0 - delta; // sorry, this is an infinite volume :(

  double equi_pdf;
  double t;
  const float u = pointsampler(p, s_dim_free_path);
  if(D > 1e-8)
  {
    double thetaA = atan(a/D);
    double thetaB = atan(b/D);

    t = D * tan((1.0 - u)*thetaA + u*thetaB);
    equi_pdf = D/fabs(thetaA - thetaB)/(D*D + t*t);
  }
  else
  {
    t = a*b / (b + (a - b)*u);
    equi_pdf = a*b / (b - a) / (t*t);
  }
  p->e[v].dist = delta + t;
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = p->e[v].vol.shader;

  // update vertex position:
  for(int k=0;k<3;k++)
    p->v[v].hit.x[k] = p->v[v-1].hit.x[k] + p->e[v].dist * p->e[v].omega[k];

  if(!path_visible(p, v  )) goto fail;

  (void)shader_prepare(p, v);

  // include distance sampling in projected solid angle pdf on vertex (will be
  // converted to on-surface on the outside, for instance in path_extend)
  p->v[v].pdf = mf_mul(p->v[v].pdf, mf_set1(equi_pdf));
  p->v[v].pdf = mf_mul(p->v[v].pdf, mf_set1(path_G(p, v)));
  p->e[v].transmittance = shader_vol_transmittance(p, v); // also inits edge transmittance and pdf, need to overwrite
  p->e[v].pdf = mf_set1(equi_pdf);

  // get the rest of the measurement contribution which we did not sample (mostly 1/mu_t):
  p->v[v].throughput = mf_mul(p->v[v].throughput, mf_div(p->e[v].transmittance, p->e[v].pdf));

  p->v[v].tech = s_tech_equiangular;
  p->length = v+1; // instruct pointsampler to get new dimensions

  // init edge v..v+1:
  p->length = v+2; // constructed vertex v and v+1
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->v[v+1].hit.x[k] - p->v[v].hit.x[k];
  p->e[v+1].dist = sqrtf(dotproduct(p->e[v+1].omega, p->e[v+1].omega));
  for(int k=0;k<3;k++) p->e[v+1].omega[k] /= p->e[v+1].dist;
  if(!path_visible(p, v+1)) goto fail;

  // we sampled this already during nee. in particular, don't use any more
  // random numbers on the last vertex (there aren't any, nee has to be
  // performed one vertex back)
  p->sensor.aperture_set = 1;
  if(p->v[0].mode & s_emit)
    edf = view_cam_connect(p);
  // TODO: also reconnect to light (in case of cosine power EDF)

  mf_t bsdf = shader_brdf(p, v); // set mode on vertex v
  if(!(p->v[v].mode & s_volume)) goto fail;
  if(mf_all(mf_eq(p->v[v].interior.mu_t, mf_set1(0.0f))))
    goto fail;

  if(path_edge_init_volume(p, v+1)) goto fail;
  if(p->e[v+1].vol.shader == -1) goto fail; // no volume no equiangular

  p->e[v+1].transmittance = shader_vol_transmittance(p, v+1);
  p->v[v+1].throughput = mf_mul(mf_mul(edf, p->e[v+1].transmittance), mf_mul(mf_set1(path_G(p, v+1)), mf_mul(p->v[v].throughput, bsdf)));

  if(0)
  {
fail:
    p->v[v].pdf = mf_set1(0.0f);
    p->throughput = mf_set1(0.0f);
    p->v[v].throughput = mf_set1(0.0f);
    p->v[v].tech = s_tech_equiangular;
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
  p->v[v+1].total_throughput = p->v[v+1].throughput; // store for path_pop()

  return 0;
}
