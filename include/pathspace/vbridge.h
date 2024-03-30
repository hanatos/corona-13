#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"

#include <math.h>

static inline int
vbridge_possible(const path_t *p, const int v, const int n)
{
  // if there are no non-specular components, fail next event estimation.
  if(!(p->v[v].material_modes & (s_diffuse | s_glossy)))
    return 0;
  if(v == 0) return 0; // as nee: do not connect light source directly
  if(v + n >= PATHSPACE_MAX_VERTS) return 0; // not enough vertices on list
  return 1;
}

// this returns the product vertex area measure pdf of sampling
// the whole bridge (including the light vertex x_n)
static inline mf_t
vbridge_pdf(
    const path_t* p, // the full path
    int v,           // the index of the light vertex x_n
    int n)           // the length of the bridge from x_0..x_n
{
  if(p->length < 3) return mf_set1(0.0f); // no next event for 2-vertex paths.
  if(n < 2) return mf_set1(0.0f); // bridge will add at least 2 vertices
  assert(v == 0 || v == p->length-1);
  int v0 = v ? v-n : v+n;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int e0 = v ? v-n+1 : n;  // towards double scattering vertex
  int en = v ? v   : 1;    // towards last point
  int vb = v ? v-n : v;    // for iteration: begin
  int ve = v ? v   : v-n;  // and end indices
  if(!vbridge_possible(p, v0, n)) return mf_set1(0.0f);
  for(int i=vb+1;i<ve;i++) // check inner vertices are in volume
    if(!primid_invalid(p->v[i].hit.prim)) return mf_set1(0.0f);

  // compute vertex area measure pdf of x_n on the light:
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

  // now compute pdf of the rest:
  // we assume fixed n, multiply P_n from outside
  // multiply all inner phase function pdf:
  float sum_d = 0.0f;
  float G = 1.0f;
  for(int i=vb+1;i<ve;i++)
  {
    res = mf_mul(res, shader_pdf(p, i));
    sum_d += p->e[i].dist;
    G *= path_G(p, i);
  }
  sum_d += p->e[ve].dist;

  const float distv[] = {
    p->v[ve].hit.x[0] - p->v[vb].hit.x[0],
    p->v[ve].hit.x[1] - p->v[vb].hit.x[1],
    p->v[ve].hit.x[2] - p->v[vb].hit.x[2]};
  uint32_t factorial = 1; // (n-1)!
  for(int i=2;i<n;i++) factorial *= i;
  const float s = sqrt(dotproduct(distv, distv));
  s = G*s*s*s * factorial / powf(sum_d, n);

  return mf_mul(res, mf_set1(s));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// next event estimation of the reverse path from v+1
static inline mf_t
vbridge_pdf_adjoint(const path_t *path, int v, int n)
{
  // luckily these guys know already about v==0 or v==length-1
  return vbridge_pdf(path, v, n);
}

// sample bridge to next event. returns != 0 on failure
static inline int
vbridge_sample(path_t *p)
{
  const int n = 2; // XXX TODO sample between 2..something
  assert(p->length >= n); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  if(p->length + n > PATHSPACE_MAX_VERTS) return 1; // need to append n new vertices
  const int v = p->length-1; // constructing new vertices starting here

  if(!vbridge_possible(p, v, n)) goto fail;

  // sample vertex vn as endpoint (lights or camera)
  int vn = v+1;
  mf_t edf = mf_set1(0.0f);
  mf_t p_egv[3];
  lights_pdf_type(p, v, p_egv);

  // constructing x_n here, will move to end of path later
  memset(p->v+vn, 0, sizeof(vertex_t));
  memset(p->e+vn, 0, sizeof(edge_t));
  // TODO: fix MIS stuff for bridges too
  // p->v[vn].pdf_mnee = mf_set1(-1.f);

  // instruct kelemen mlt to use new random numbers:
  p->v[vn].rand_beg = p->v[v].rand_beg + p->v[v].rand_cnt;

  // set technique
  p->v[vn].tech = s_tech_vbridge;

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
    p->v[vn].pdf = mf_mul(p->v[vn].pdf, p_egv[0]);
  }
  else if(rand < mf(p_egv[0],0)+mf(p_egv[1],0))
  { // connect to geo lights
    edf = lights_sample_next_event(p);
    // compensate for envmap sampling probability
    p->v[vn].pdf = mf_mul(p->v[vn].pdf, p_egv[1]);
    edf = mf_div(edf, mf_set1(mf(p_egv[1], 0)));
  }
  else if(mf(p_egv[2],0) > 0)
  { // connect to volume lights
    edf = light_volume_sample_nee(p, v);
    if(!mf_any(mf_gt(edf, mf_set1(.0f)))) goto fail;
    p->v[vn].pdf = mf_mul(p->v[vn].pdf, p_egv[2]);
    edf = mf_div(edf, mf_set1(mf(p_egv[2],0)));
    // init segment
    for(int k=0;k<3;k++)
      p->e[vn].omega[k] = p->v[vn].hit.x[k] - p->v[v].hit.x[k];
    p->e[vn].dist = sqrtf(dotproduct(p->e[vn].omega, p->e[vn].omega));
    for(int k=0;k<3;k++) p->e[vn].omega[k] *= 1./p->e[vn].dist;
  }

  // ask edf and bsdf for their consent:
  if(!mf_any(mf_gt(edf, mf_set1(0.0f)))) goto fail;
  // compute brdf and throughput (will ignore, just check for specular)
  mf_t bsdf = shader_brdf(p, v); // also set mode on vertex v
  if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.

  // determine side of surface and volume from that (brdf sets mode)
  if(path_edge_init_volume(p, vn)) goto fail;

  // move vertex vn to end of list, free up space for x_1..x_{n-1}, our in-betweens
  p->e[v+n] = p->e[vn];
  p->v[v+n] = p->v[vn];
  p->v[v+n].rand_beg = p->v[v].rand_beg + s_dim_num_nee;
  vn = v+n; // now the index is pointing to the vertex at the correct position

  // get volume properties:
  if(p->e[vn].vol.shader == -1) goto fail; // no volume

  // w_1 on e[v+1] is now inited to point directly to the nee vertex.
  // vertex_t xn = p->v[vn]; // store target away, we will overwrite it
  const float to_target_xn[] = {
    p->v[vn].hit.x[0] - p->v[v].hit.x[0],
    p->v[vn].hit.x[1] - p->v[v].hit.x[1],
    p->v[vn].hit.x[2] - p->v[v].hit.x[2]};

  // sample distances and phase functions
  for(int i=v+1;i<=n;i++)
  {
    p->length = i; // instruct pointsampler to get new dimensions
    // TODO: use homogeneous medium scattering code here! use constant parameters
    shader_vol_sample(p, i); // writes e[i].dist
    for(int k=0;k<3;k++) // update hit position:
      p->v[i].hit.x[k] = p->v[i-1].x[k] + p->e[i].dist * p->e[i].omega[k];
    p->v[i].rand_beg = p->v[i+1].rand_beg + s_dim_num_extend;
    shader_prepare(p, i);
    p->length = i+1;
    if(i < n)
    {
      shader_sample(p);        // writes p->e[i+1].omega
      if(path_edge_init_volume(p, i+1)) goto fail;
      if(p->e[i+1].vol.shader == -1)    goto fail; // no volume
    }
  }
  
  // setup quaternion rotation
  const float to_trace_xn[] = {
    p->v[vn].hit.x[0] - p->v[v].hit.x[0],
    p->v[vn].hit.x[1] - p->v[v].hit.x[1],
    p->v[vn].hit.x[2] - p->v[v].hit.x[2]};
  const float len_target = sqrtf(dotproduct(to_target_xn,to_target_xn));
  const float len_traced = sqrtf(dotproduct(to_trace_xn,to_trace_xn));
  for(int k=0;k<3;k++) to_trace_xn[k] /= len_traced;
  for(int k=0;k<3;k++) to_target_xn[k] /= len_target;
  quaternion_t q;
  crossproduct(to_trace_xn, to_target_xn, q.x);
  q.w = 1.0f + dotproduct(to_trace_xn, to_target_xn);
  const float l = sqrtf(q.x[0]*q.x[0] + q.x[1]*q.x[1] + q.x[2]*q.x[2] + q.w*q.w);
  for(int k=0;k<3;k++) q.x[k] /= l;
  q.w /= l;

  // rotate/scale vertices
  float sum_d = 0.0f;
  for(int i=v+1;i<=n;i++)
  {
    // rotate vertex
    float x[] = {
      p->v[i].hit.x[0] - p->v[v].hit.x[0],
      p->v[i].hit.x[1] - p->v[v].hit.x[1],
      p->v[i].hit.x[2] - p->v[v].hit.x[2]};
    quaternion_transform(&q, &x);
    for(int k=0;k<3;k++) p->v[i].hit.x[k] = 
      p->v[v].hit.x[k] + len_target/len_traced * x[k];
    // adjust edge
    quaternion_transform(&q, p->e[i].omega);
    p->e[i].dist *= len_target/len_traced;
    p->v[i].pdf = 1.0f; // will include it all in the last vertex
    sum_d += p->e[i].dist;
  }

  // compute throughput and pdf
  mf_t f = path_measurement_contribution_dwp(p, v, vn);
  mf_t f = shader_brdf(p, v); // start vertex
  for(int i=v+1;i<n;i++)
  {
    f = mf_mul(f, shader_vol_transmittance(p, i));
    f = mf_mul(f, path->v[i].vol.mu_s);// XXX ???
  }
  f = mf_mul(f, shader_vol_transmittance(p, nv));
  f = mf_mul(f, lights_eval_vertex(p, vn)); // end vertex
  uint32_t factorial = 1; // (n-1)!
  for(int i=2;i<n;i++) factorial *= i;
  const float s = len_target;
  s = s*s*s * factorial / powf(sum_d, n);
  mf_t pdf = mf_set1(s);
  p->v[vn].throughput = mf_div(f, pdf);
  p->v[vn].pdf = vbridge_pdf(p, vn, n); // area measure

  if(0)
  {
fail:
    p->v[vn].pdf = mf_set1(0.0f);
    p->throughput = mf_set1(0.0f);
    p->v[vn].throughput = mf_set1(0.0f);
    p->v[vn].flags = s_none;
    p->v[vn].mode = s_absorb;
    p->v[vn].rand_cnt = s_dim_num_nee;
    p->v[vn] = p->v[v];
    // add enough to be popped off
    p->length = v+n+1; // constructed vertices (even if they absorb)
    return 0;
  }
  p->v[v].rand_cnt = s_dim_num_nee;
  for(int i=1;i<=n;i++)
    p->v[v+i].rand_cnt = s_dim_num_extend;
  p->length = v+n+1;

  p->v[vn].total_throughput = p->v[vn].throughput; // store for path_pop()
  return 0;
}
