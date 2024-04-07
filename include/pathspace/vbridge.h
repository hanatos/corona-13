#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"
#include "points.h"

#include <math.h>

static mf_t num_verts_P(const float dist, const mf_t mu_t, int n);
static inline int 
num_verts_sample(const float dist, const mf_t mu_t, mf_t *P)
{
#if 0
  if(P) P[0] = 1.0f; // deterministic
  return 3;// XXX
  float meand = 1.0f/HOMO_MU_T;
  return (int)(dist / meand);
#else
  // poisson
  float T = 0.0f;
  int n = -1;
  const float l = mf(mu_t, 0) * dist;
  while(T <= 1.0f)
  {
    float xi = points_rand(rt.points, common_get_threadid());
    T -= log(1.0f-xi)/l;
    n++;
  }
  if(P) P[0] = num_verts_P(dist, mu_t, n);
  return n;
#endif
}

static inline mf_t
num_verts_P(const float dist, const mf_t mu_t, int n)
{
#if 0
  if(num_verts_sample(dist, 0) != n) return 0.0f;
  return 1.0f; // deterministic
#else
  const mf_t l = mf_mul(mu_t, mf_set1(dist));
  uint32_t fac = 1; // n!
  for(int i=2;i<=n;i++) fac *= i;
  float res[MF_COUNT];
  // sucks, but clang is unhappy about non-const indexing in a loop
  res[0] = powf(mf(l, 0), n) * expf(-mf(l, 0)) / fac;
#if MF_COUNT > 1
  res[1] = powf(mf(l, 1), n) * expf(-mf(l, 1)) / fac;
  res[2] = powf(mf(l, 2), n) * expf(-mf(l, 2)) / fac;
  res[3] = powf(mf(l, 3), n) * expf(-mf(l, 3)) / fac;
#if MF_COUNT > 4
  res[4] = powf(mf(l, 4), n) * expf(-mf(l, 4)) / fac;
  res[5] = powf(mf(l, 5), n) * expf(-mf(l, 5)) / fac;
  res[6] = powf(mf(l, 6), n) * expf(-mf(l, 6)) / fac;
  res[7] = powf(mf(l, 7), n) * expf(-mf(l, 7)) / fac;
#endif
#endif
  return mf_loadu(res);
#endif
}

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
  if(p->length < 4) return mf_set1(0.0f); // no next event for 2-vertex paths.
  if(n < 2) return mf_set1(0.0f); // bridge will add at least 2 vertices
  // assert(v == 0 || v == p->length-1);
  if(!(v == 0 || v == p->length-1)) return mf_set1(0.0f);
  int v0 = v ? v-n : v+n;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int vb = v ? v-n : v;    // for iteration: begin
  int ve = v ? v   : v-n;  // and end indices
  if(!vbridge_possible(p, v0, n)) return mf_set1(0.0f);
  for(int i=vb+1;i<ve;i++) // check inner vertices are in volume
    if(!primid_invalid(p->v[i].hit.prim)) return mf_set1(0.0f);

  // compute vertex area measure pdf of x_n on the light:
  mf_t p_egv[3];
  lights_pdf_type(p, v0, p_egv);

  mf_t res = mf_set1(0.0f);
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
  // n is fixed, multiply P_n from outside
  // multiply all inner phase function pdf:
  double sum_d = 0.0;
  double G = 1.0;
  for(int i=vb+1;i<=ve;i++)
  {
    if(i<ve) res = mf_mul(res, shader_pdf(p, i));
    sum_d += p->e[i].dist;
    G *= path_G(p, i);
  }
  sum_d += p->e[ve].dist;

  const float distv[] = {
    p->v[ve].hit.x[0] - p->v[vb].hit.x[0],
    p->v[ve].hit.x[1] - p->v[vb].hit.x[1],
    p->v[ve].hit.x[2] - p->v[vb].hit.x[2]};
  double factorial = 1; // (n-1)!
  for(int i=2;i<n;i++) factorial *= i;
  double s = sqrt(dotproduct(distv, distv));
  const mf_t P_n = num_verts_P(s, p->v[vb].interior.mu_t, n);
  s = G*s*s*s * factorial / pow(sum_d, n);

  return mf_mul(res, mf_mul(P_n, mf_set1(s)));
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
  assert(p->length >= 2); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  const int v = p->length-1; // constructing new vertices starting here

  // sample vertex vn as endpoint (lights or camera)
  int vn = v+1;
  mf_t edf = mf_set1(0.0f);
  mf_t p_egv[3];
  lights_pdf_type(p, v, p_egv);

  // constructing x_n here, will move to end of path later
  memset(p->v+vn, 0, sizeof(vertex_t));
  memset(p->e+vn, 0, sizeof(edge_t));

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
    p->v[vn].pdf = mf_mul(p->v[vn].pdf, p_egv[2]);
    edf = mf_div(edf, mf_set1(mf(p_egv[2],0)));
    // init segment
    for(int k=0;k<3;k++)
      p->e[vn].omega[k] = p->v[vn].hit.x[k] - p->v[v].hit.x[k];
    p->e[vn].dist = sqrtf(dotproduct(p->e[vn].omega, p->e[vn].omega));
    for(int k=0;k<3;k++) p->e[vn].omega[k] *= 1./p->e[vn].dist;
  }

  // w_1 on e[v+1] is now inited to point directly to the nee vertex.
  float to_target_xn[] = {
    p->v[vn].hit.x[0] - p->v[v].hit.x[0],
    p->v[vn].hit.x[1] - p->v[v].hit.x[1],
    p->v[vn].hit.x[2] - p->v[v].hit.x[2]};
  const float len_target = sqrtf(dotproduct(to_target_xn,to_target_xn));

  mf_t P_n = mf_set1(0.0f);
  const int n = num_verts_sample(len_target, p->v[v].interior.mu_t, &P_n);
  // fprintf(stderr, " v %d n %d vn %d\n", v, n, vn);
  if(n < 2) return 1;
  if(p->length + n > PATHSPACE_MAX_VERTS) return 1; // need to append n new vertices

  // fprintf(stderr, "edf %g\n", mf(edf, 0));
  // ask edf and bsdf for their consent:
  // if(!mf_any(mf_gt(edf, mf_set1(0.0f)))) goto fail; // this is the wrong direction
  // compute brdf and throughput (will ignore, just check for specular)
  mf_t bsdf = shader_brdf(p, v); // also set mode on vertex v
  (void)bsdf;
  // fprintf(stderr, "brdf %g\n", mf(bsdf, 0));
  // if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials. // need to check flags instead

  // fprintf(stderr, "XXX v %d n %d vn %d\n", v, n, vn);
  // determine side of surface and volume from that (brdf sets mode)
  if(path_edge_init_volume(p, vn)) goto fail;
  // fprintf(stderr, "III v %d n %d vn %d\n", v, n, vn);

  // p->length = v+2;
  // path_print(p, stderr);
  // move vertex vn to end of list, free up space for x_1..x_{n-1}, our in-betweens
  p->e[v+n] = p->e[vn];
  p->v[v+n] = p->v[vn];
  memset(p->v+vn, 0, sizeof(vertex_t));
  p->v[v+n].rand_beg = p->v[v].rand_beg + s_dim_num_nee;
  vn = v+n; // now the index is pointing to the vertex at the correct position

  // get volume properties:
  if(p->e[vn].vol.shader == -1) goto fail; // no volume

  // sample distances and phase functions
  for(int i=v+1;i<=vn;i++)
  {
    p->length = i; // instruct pointsampler to get new dimensions
    p->e[i].dist = -logf(1.0f-pointsampler(p, s_dim_free_path));// mu_t cancels anyways
    for(int k=0;k<3;k++) // update hit position:
      p->v[i].hit.x[k] = p->v[i-1].hit.x[k] + p->e[i].dist * p->e[i].omega[k];
    if(i < vn)
    {
      p->v[i].hit.prim = INVALID_PRIMID;
      p->v[i].hit.shader = p->e[i].vol.shader;
    }
    p->v[i].rand_beg = p->v[i-1].rand_beg + s_dim_num_extend;
    shader_prepare(p, i);
    p->length = i+1;
    if(i < vn)
    {
      shader_sample(p);        // writes p->e[i+1].omega and sets v[i].mode
      if(path_edge_init_volume(p, i+1)) goto fail;
      if(p->e[i+1].vol.shader == -1)    goto fail; // no volume
    }
    // fprintf(stderr, "v %d mode %d vol shader %d\n", i, p->v[i].mode, p->v[i].interior.shader);
  }
  // path_print(p, stderr);
  
  // setup quaternion rotation
  float to_trace_xn[] = {
    p->v[vn].hit.x[0] - p->v[v].hit.x[0],
    p->v[vn].hit.x[1] - p->v[v].hit.x[1],
    p->v[vn].hit.x[2] - p->v[v].hit.x[2]};
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
  for(int i=v+1;i<=vn;i++)
  {
    // rotate vertex
    float x[] = {
      p->v[i].hit.x[0] - p->v[v].hit.x[0],
      p->v[i].hit.x[1] - p->v[v].hit.x[1],
      p->v[i].hit.x[2] - p->v[v].hit.x[2]};
    quaternion_transform(&q, x);
    for(int k=0;k<3;k++) p->v[i].hit.x[k] = 
      p->v[v].hit.x[k] + len_target/len_traced * x[k];
    // adjust edge
    quaternion_transform(&q, p->e[i].omega);
    p->e[i].dist *= len_target/len_traced;
    p->v[i].pdf = mf_set1(1.0f); // will include it all in the last vertex
    sum_d += p->e[i].dist;
    if(!path_visible(p, i)) goto fail;
  }

  // compute throughput and pdf
  // uint32_t factorial = 1; // (n-1)!
  // for(int i=2;i<n;i++) factorial *= i;
  // hope for better precision if we don't blow number range up so much:
  double fac = 1.0/(sum_d * sum_d);
  for(int i=2;i<n;i++) fac *= i/sum_d;
  double s = len_target;
  s = s*s*s * fac; // factorial / powf(sum_d, n);
  // fprintf(stderr, "pdf len %g, P_n %g, sum_d %g, fac %g\n", len_target, P_n, sum_d, fac);
  md_t pdf = md_mul(mf_2d(P_n), md_set1(s));

  md_t f = mf_2d(shader_brdf(p, v)); // start vertex
  for(int i=v+1;i<v+n;i++)
  {
  // fprintf(stderr, "i %d len %d, f %g, pdf %g\n", i, vn+1, md(f, 0), md(pdf, 0));
    f = md_mul(f, mf_2d(shader_vol_transmittance(p, i)));
    f = md_mul(f, mf_2d(p->v[i].interior.mu_s));
  }
  f = md_mul(f, mf_2d(shader_vol_transmittance(p, vn)));
  f = md_mul(f, mf_2d(lights_eval_vertex(p, vn))); // end vertex
  // f = md_set1(1);
  // pdf = md_set1(1);
  p->v[vn].throughput = md_2f(md_div(f, pdf));
  p->length = vn+1;
  p->v[vn].pdf = vbridge_pdf(p, vn, n); // area measure
  p->v[v].rand_cnt = s_dim_num_nee;
  for(int i=1;i<=n;i++) p->v[v+i].rand_cnt = s_dim_num_extend;

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
    p->throughput = mf_set1(0.0);
    return 1;
  }
  p->v[v].total_throughput = p->v[v].throughput;
  p->throughput = p->v[v].throughput * p->v[vn].throughput;
  return 0;
}

static inline void
vbridge_pop(path_t *p, int old_length)
{
  p->throughput = mf_set1(0.0);
  p->v[old_length].throughput = p->v[old_length].total_throughput;
  p->length = old_length;
}
