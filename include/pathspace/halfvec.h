#ifndef _PATHSPACE_HALFVEC_H
#define _PATHSPACE_HALFVEC_H

// some halfvector utility functions which we can make good use of even outside MLT

#include "pathspace.h"
#include "shader.h"
#include "lights.h"
#include "accel.h"
#include "points.h"
#include "pathspace/manifold.h"
#include "spectrum.h"
#include "pathspace/raydifferentials.h"
#include "matrix2.h"
#include "view.h"

// step how many pixels (in percent of the image width):
#define HALFVEC_MUTATION_STEP (2.f*MIN(view_width(),view_height())/100.0f)
// adjust optimal bsdf step by this amount in beckmann space:
#define HALFVEC_BSDF_STEP 1.0f
// halfvec distance perturbation step size in multiples of mean free path:
#define HALFVEC_DIST_STEP 0.1f
#define HALFVEC_BECKMANN_MIN 1e-8
#define HALFVEC_BECKMANN_MAX 1.7
#define HALFVEC_REL_SPATIAL_EPS 1e-3f
#define HALFVEC_SQR_HV_EPS 1e-8f


// DEBUG
#define USE_RAYDIFF 1

typedef struct halfvec_stats_t
{
  uint64_t project_failed;
  uint64_t propagate_failed;
  uint64_t perturb_failed;
  uint64_t no_throughput;
  uint64_t mutations;
  uint64_t newton_walks;
  uint64_t successful_newton_walks;
  uint64_t mutations_proposed;
  uint64_t reverse_check_failed;
  uint64_t sum_iterations;
  uint64_t sum_successful_reduces;
  uint64_t raydiff_compute_failed;
  uint64_t error_increased;
}
halfvec_stats_t;

static inline float halfvec_to_worldspace(
    halfvec_stats_t *d,  // stats struct
    path_t *tent,        // tentative path to be updated
    const float *h,      // e-s+1 target half vectors and distance constraints (including a slot for v[s] and v[e])
    const int s,         // fixed vertex at start of chain
    const int e)         // fixed vertex at end of chain
{
  d->newton_walks++;
  const int n = e-s;
  float dh[(n+1)*3];
  float err = FLT_MAX;
  const int num_it = 15; // maximum number of iterations before giving up
  float newerr = 0.0f;
  int reduce_it = 0;
  int err_inc = 0;
  float step = 1.0f;
  for(int k=0;k<3;k++) dh[k] = 0;
  for(int k=0;k<3;k++) dh[3*n+k] = 0;
  for(int i=1;i<n;i++)
  {
    dh[3*i+0] = h[3*i+0] - tent->v[s+i].diffgeo.h[0];
    dh[3*i+1] = h[3*i+1] - tent->v[s+i].diffgeo.h[1];
    if(tent->v[s+i].mode & s_volume) dh[3*i+2] = h[3*i+2] - tent->v[s+i].diffgeo.h[2];
    else dh[3*i+2] = 0.0f;
  }
  path_t backup = *tent;
  for(int it=0;it<num_it;it++)
  {
    // update vertex positions using first-order approximation:
    if(manifold_map_h_to_x(tent, tent, dh, step, s, e)) return 0.0f;

    // trace tent path towards all updated vertex positions:
    for(int v=s+1;v<=e;v++)
    {
      tent->e[v].transmittance = mf_set1(0.0f);
      if(v != e || !(tent->v[e].flags & s_environment))
      { // keep projected position at v[v]
        if(path_project(tent, v, s_propagate_mutate)) goto reduce_stepsize;
      }
      else
      { // or the direction to the envmap
        if(path_propagate(tent, v, s_propagate_mutate)) goto reduce_stepsize;
      }
      // check integrity of the new vertex
      if((tent->v[v-1].mode != backup.v[v-1].mode) ||
         (tent->v[v].flags  != backup.v[v].flags) ||
         (tent->v[v].hit.shader != backup.v[v].hit.shader) ||
         (tent->v[v].interior.shader != backup.v[v].interior.shader) ||
         (primid_invalid(tent->v[v].hit.prim) != primid_invalid(backup.v[v].hit.prim)))
        goto reduce_stepsize;
      // check integrity of the sidedness on the last vertex:
      if((tent->v[v-1].mode & (s_reflect | s_transmit)) && v > 1 &&
          ((dotproduct( tent->e[v].omega,  tent->v[v-1].hit.n) > 0.0f) !=
           (dotproduct(backup.e[v].omega, backup.v[v-1].hit.n) > 0.0f)))
        goto reduce_stepsize;
    }

    // re-initialize transport matrices to match new half vectors and positions
    const double det = manifold_compute_tangents(tent, s, e);
    if(det == 0.0) return 0.0f;
    // compute residual relative half vector distance
    newerr = 0.0f;
    for(int i=1;i<n;i++)
    {
      dh[3*i+0] = h[3*i+0] - tent->v[s+i].diffgeo.h[0];
      dh[3*i+1] = h[3*i+1] - tent->v[s+i].diffgeo.h[1];
      if(tent->v[s+i].mode & s_volume) dh[3*i+2] = h[3*i+2] - tent->v[s+i].diffgeo.h[2];
      else dh[3*i+2] = 0.0f;
      // compute error on projected half vector disk:
      const float c2 = 1.0f/sqrtf(1.0f + h[3*i]*h[3*i] + h[3*i+1]*h[3*i+1]);
      const float c0 = h[3*i]*c2, c1 = h[3*i+1]*c2;
      const float h2 = 1.0f/sqrtf(1.0f + tent->v[s+i].diffgeo.h[0]*tent->v[s+i].diffgeo.h[0] + tent->v[s+i].diffgeo.h[1]*tent->v[s+i].diffgeo.h[1]);
      const float h0 = tent->v[s+i].diffgeo.h[0]*h2, h1 = tent->v[s+i].diffgeo.h[1]*h2;
      // put error close to [0,1] range by normalising out the target value.
      const float dht = tent->v[s+i].mode & s_volume ? (tent->v[s+i].diffgeo.h[2] - h[3*i+2])/fabsf(h[3*i+2]) : 0.0f;
      newerr = fmaxf(newerr, (h0-c0)*(h0-c0) + (h1-c1)*(h1-c1) + dht*dht);
    }

    // converged?
    if(newerr < HALFVEC_SQR_HV_EPS)
    {
      d->sum_iterations += it+1;
      d->successful_newton_walks ++;
      d->sum_successful_reduces += reduce_it;
      d->error_increased += err_inc;
      return 1.0f;
    }
    else if(newerr > err)
    {
      err_inc ++;
reduce_stepsize:
      reduce_it ++;
      *tent = backup;
      step *= 0.5f;
      if(step < 1e-5f) break; // counts as max iteration count hit.
    }
    else
    {
      backup = *tent;
      err = newerr;
      step = fminf(1.0f, step * 2.0f);
    }
  }
  // max iterations reached, reduced stepsize sometimes
  return 0.0f;
}

// helper for halfvec_measurement. it computes |do^perp/dh^par| for each inner vertex.
static inline double _halfvec_measurement_vertex(const path_t *path, int v)
{
  double f = 1.0;
  // specular surface reflection returns its bsdf in projected half vector space already,
  // (i.e. only fresnel etc) so all the jacobians cancel out.
  if(!(path->v[v].mode & s_specular))
  {
    // cosine ratio to pretend we're computing the same measurement as the path
    // tracer (incoming geo normal, outgoing shading normal)
    assert(path->v[0].mode & s_sensor); // only works for pt direction in this writing:
    f = MIN(4.0f, fabs(dotproduct(path->v[v].hit.gn, path->e[v].omega) /
        MAX(1e-8f, fabs(dotproduct(path->v[v].hit.n, path->e[v].omega)))));
    // |dw^perp_k/dh^par_k|
    float h[3], H[3]; // reconstruct world space half vector
    // plane/plane half vector lives in ortho normal basis hit.{a,b,n}
    // (and usually points into hemisphere of optically thinner medium, but we don't care)
    h[0] = path->v[v].diffgeo.h[0];
    h[1] = path->v[v].diffgeo.h[1];
    h[2] = 1.0f/sqrtf(1.0 + h[0]*h[0] + h[1]*h[1]);
    h[0] *= h[2];
    h[1] *= h[2];
    // h points into random hemisphere, but sign of H
    // doesn't matter further down in the computation
    for(int k=0;k<3;k++)
      H[k] = path->v[v].hit.a[k] * h[0] +
        path->v[v].hit.b[k] * h[1] +
        path->v[v].hit.n[k] * h[2];

    if(path->v[v].mode & s_reflect)
    {
      f *= dotproduct(path->v[v].hit.n, path->e[v+1].omega) * h[2]*h[2]*h[2];
      f *= 4.0f * dotproduct(H, path->e[v+1].omega);
    }
    else if(path->v[v].mode & s_transmit)
    {
      float eta = mf(path_eta_ratio(path, v), 0);
      if(eta != 1.0f)
      {
        f *= dotproduct(path->v[v].hit.n, path->e[v+1].omega) * h[2]*h[2]*h[2];
        eta = 1./eta;
        const float dot_H_wo = dotproduct(H, path->e[v+1].omega);
        if(fabsf(dot_H_wo) < 1e-4) return 0.0;
        float den = -dotproduct(H, path->e[v].omega) +
          eta*dotproduct(H, path->e[v+1].omega);
        f *= den*den/(eta*eta*dot_H_wo);
      }
      // else do/dh == 1
    }
    else if((path->v[v].mode & s_volume) || (path->v[v].mode & s_fiber))
    { // verified in volumes against f(X)/J:
      float h[3] = {path->v[v].diffgeo.h[0], path->v[v].diffgeo.h[1], 1.0};
      const float z = 2.0f/(h[0]*h[0] + h[1]*h[1] + 1.0f);
      // no *=, volumes don't need the cosine factor.
      f = z*z; // half vector now in 2d plane after transform via riemannian sphere
      // TODO: fibers probably need a sin(tangent, omega) here. note
      // that this code should not currently be called on fibers, as
      // we don't have a meaningful tangent frame derivative either.
    }
  }
  return fabs(f);
}

static inline md_t _halfvec_measurement(path_t *path, const int s, const int e)
{
#define PER_VERTEX_JACOBIAN f = md_mul(f, md_set1(_halfvec_measurement_vertex(path, v)));
#define MEASUREMENT_BEG s
#define MEASUREMENT_END e
#include "pathspace/measurement.h"
}

// compute measurement contribution in half vector space (i.e. f * |J| if f is in vertex area measure)
// compute tangents has to be called already before calling this.
// essentially a faster and more precise version of:
// return path_measurement_contribution_dx(path) / path->cache.dh_dx;
static inline md_t halfvec_measurement(path_t *path, const int s, const int e)
{
  // compute generalized geometric term for full path
  // just use transfer matrix determinant directly. dpdu/dpdv are orthonormal, so
  // mapping them is just wasted compute cycles.
  float dw0_dxn = path_G(path, s+1); // = dw0_dx1;
  // if(e-s > 1) dw0_dxn *= mat2_det(path->v[s+1].diffgeo.Tp); // *= dx1_dxn;
  if(e-s > 1)
  {
    // const float dx1_dxn = fminf(4.0, fabsf(mat2_det(path->v[s+1].diffgeo.Tp)));
    const float dx1_dxn = mat2_det(path->v[s+1].diffgeo.Tp);
    dw0_dxn *= dx1_dxn;
    // account for compression of ray density due to refractive indices.
    // we ignore this everywhere and that actually leads to correct results. hooray :)
    // const float n1 = path->e[s+1].vol.ior;
    // const float n2 = path->e[e].vol.ior;
    // const float f = (n1/n2)*(n1/n2);
    // dw0_dxn *= f;
  }
  // TODO: do we need mf_abs?
  return md_mul(md_set1(fabs(dw0_dxn)), _halfvec_measurement(path, s, e));
}

static inline float _halfvec_dist_stepsize(
    const path_t *path,
    const int v)
{
  if(path->v[v].mode & s_volume)
    return HALFVEC_DIST_STEP/mf(path->v[v].interior.mu_t, 0);
  return 1.0f;
}

// returns resulting step size in ray differential space (beckmann + anisotropic pixel step vectors)
static inline float _halfvec_bsdf_stepsize(
    const path_t *path,    // current path
    const int v)           // vertex number in question
{
  if(path->v[v].mode & s_volume)
  {
    // const float *h = path->v[v].diffgeo.h;
    // const float z = 2.0f/(h[0]*h[0] + h[1]*h[1] + 1.0f);

    // convert mean cosine to mean sine and that to riemannian plane
    // and convert that to standard deviation using the same reasoning
    // as a MAD estimate (*1.4826):
    const float g = path->v[v].interior.mean_cos;
    // return HALFVEC_BSDF_STEP * sqrtf(1.0f - g*g) * 1.4826f / z;
    // (seems to be way too big, try something smaller:)
    // return HALFVEC_BSDF_STEP * sqrtf(1.0f - g*g) * 0.001f / z;
    return HALFVEC_BSDF_STEP * sqrtf(1.0f - g*g) * 0.001f;
  }
  // optimal stepsize for beckmann space bsdf bandwidth:
  return CLAMP(HALFVEC_BSDF_STEP * path->v[v].shading.roughness * sqrtf(2.0f/M_PI), HALFVEC_BECKMANN_MIN, HALFVEC_BECKMANN_MAX);
}

static inline void _halfvec_compute_stepsizes(const path_t *path, float *R, float *rd_u, float *rd_v, const float *rd_i, const float *rd_j, const int n)
{
  // fill trivial matrices and eigenvalues for inner vertices outside the sensor connecting segment:
  float *Rk = R + 9*n;
  for(int k=n;k<path->length-1;k++)
  {
    mat3_set_identity(Rk);
    rd_v[k] = rd_u[k] = 1.0; // means use bsdf/distance proposal directly
    Rk += 9;
  }

  // distribute ray differential step sizes to (stochastically) sum up to one pixel step.
  // determine ratio by roughness
  // XXX TODO: check whether mu_t / distance step can participate in this!
  float stepsizes[PATHSPACE_MAX_VERTS] = {0.0f};
  float sum = 0.f;
  for(int k=1;k<n;k++)
    sum += (stepsizes[k] = _halfvec_bsdf_stepsize(path, k));
  for(int k=1;k<n;k++)
    stepsizes[k] /= sum;

  Rk = R + 9;
  for(int k=1;k<n;k++)
  {
#if USE_RAYDIFF != 1
    // DEBUG: switch off ray diffs:
    goto no_raydiff;
#endif
    // transform offset at x[1] for a step in pixel i and pixel j coordinates to halfvector at vertex k:
    float hu[3], hv[3], ht[3];
    const float rd_k[] = {0.0f, 0.0f, 1.0f};
    mat3_mulv(Rk, rd_i, hu);
    mat3_mulv(Rk, rd_j, hv);
    mat3_mulv(Rk, rd_k, ht);

    if(!(hu[0] == hu[0])|| !(hu[1] == hu[1])|| !(hu[2] == hu[2])||
       !(hv[0] == hv[0])|| !(hv[1] == hv[1])|| !(hv[2] == hv[2])||
       !(ht[0] == ht[0])|| !(ht[1] == ht[1])|| !(ht[2] == ht[2])) goto no_raydiff;
#if 0
      mat3_set_identity(Rk);
    for(int i=0;i<3;i++) Rk[0+i] = hu[i];
    for(int i=0;i<3;i++) Rk[3+i] = hv[i];
    for(int i=0;i<3;i++) Rk[6+i] = ht[i];
      rd_u[k] = 1.0f;
      rd_v[k] = 1.0f;
#else

    // now we're actually interested in the ellipse in the 2d plane through our
    // 3d parameter space that corresponds to the pixel footprint.

    // transform into bsdf_u bsdf_v distance_w space (anisotropic scale) and clip here.
    // S scales the half vector 3d cube to unit variance along all constrant axes:
    //     | 1/sbu            |
    // S = |       1/sbv      |
    //     |             1/sd |
    const float v_bsdf0 = _halfvec_bsdf_stepsize(path, k);
    const float v_bsdf1 = v_bsdf0; // TODO: support anisotropic bsdf via two roughnesses here
    const float v_dist = _halfvec_dist_stepsize(path, k);
    const float s[3] = {v_bsdf0, v_bsdf1, v_dist};
    for(int k=0;k<3;k++) hu[k] /= s[k];
    for(int k=0;k<3;k++) hv[k] /= s[k];
    for(int k=0;k<3;k++) ht[k] /= s[k];

    float hn[3], a[3], b[3];
    // hn == ht and ortho normalise the other two before svd
    // crossproduct(hu, hv, hn);
    for(int k=0;k<3;k++) hn[k] = ht[k];
    normalise(hn);
    if(!(hn[0] == hn[0])||
       !(hn[1] == hn[1])||
       !(hn[2] == hn[2])) goto no_raydiff;
    // for(int k=0;k<3;k++) a[k] = hu[k];
    // normalise(a);
    // crossproduct(hn, a, b);
    get_onb(hn, a, b);
    // R3 rotates normalised hn = (0,0,1), R2 leaves third dimension alone:
    //       | a b hn |      | . . 0 |
    //  R3 = | a b hn | R2 = | . . 0 |  => R3 * R2 transforms single pixel steps in screen space to constraint offsets at vertex v[k].
    //       | a b hn |      | 0 0 1 |

    // construct R2 by projecting hu, hv into a, b basis:
    Rk[0] = dotproduct(hu, a);
    Rk[2] = dotproduct(hu, b);
    Rk[1] = dotproduct(hv, a);
    Rk[3] = dotproduct(hv, b);

    // the transformed basis vectors (column vectors
    // of R_k) are non-orthogonal and not normalised. plus, we'd like to scale them
    // isotropically by a certain percentage of the image width for steerable stratification
    // (this is the 2% parameter). this scaling will only affect overall size, not anisotropy:
    //
    const float iso_scale = HALFVEC_MUTATION_STEP;
    //
    // for the anisotropy, we'd like to know the main axes of the ellipse (one std dev from
    // the center), because we can adjust this later on according to the bsdf without
    // accidentally rotating the main axes again.
    //
    // perform an incomplete 2x2 svd on Rk:
    // find U such that U * Rk * Rk' * U'=diag
    // Su = Rk*Rk'
    float Rkp[4], Su[4];
    mat2_transpose(Rk, Rkp);
    mat2_mul(Rk, Rkp, Su);

    // this works, but involves cumbersome trig:
    const float phi = - .5f * atan2f(Su[1] + Su[2], Su[0] - Su[3]);
    float cos_phi, sin_phi;
    sincosf(phi, &sin_phi, &cos_phi);

    // write back new rotation matrix (= U')
    // R2 = |  cos_phi sin_phi |
    //      | -sin_phi cos_phi |
    for(int i=0;i<3;i++) Rk[0+i] =  a[i] * cos_phi + b[i] * sin_phi;
    for(int i=0;i<3;i++) Rk[3+i] = -a[i] * sin_phi + b[i] * cos_phi;
    for(int i=0;i<3;i++) Rk[6+i] = hn[i];

    // find the singular values from U
    const float Su_sum = Su[0] + Su[3];
    const float Su_dif = sqrtf(fmaxf(0.0, (Su[0]-Su[3])*(Su[0]-Su[3]) + 4.0f*Su[1]*Su[2]));
    // now the two axes will have anisotropic step lengths.
    // use the minimum of local distance/bsdf steps (unit sphere in this space) and ray differentials:
    rd_u[k] = MIN(1.0f, stepsizes[k] * iso_scale * sqrtf(MAX(1e-10f, (Su_sum + Su_dif)*.5f)));
    rd_v[k] = MIN(1.0f, stepsizes[k] * iso_scale * sqrtf(MAX(1e-10f, (Su_sum - Su_dif)*.5f)));
    // and rd_w[k] = 1.0f;

    assert(rd_u[k] == rd_u[k]);
    assert(rd_v[k] == rd_v[k]);
#endif
    if(0)
    {
no_raydiff:
      mat3_set_identity(Rk);
      rd_u[k] = 1.0f;
      rd_v[k] = 1.0f;
    }
    Rk += 9;
  }
}

// return dh_dx
static inline double halfvec_compute_raydifferentials(
    halfvec_stats_t *stats,
    path_t *p,              // path to compute rd for
    const int s,            // fixed vertex at start of specular chain
    const int e,            // fixed vertex at end of specular chain
    float *R,               // store ray diffs or 0
    float *rd_u,
    float *rd_v)
{
  // get the x[1] -> h[i] matrix for every vertex i:
  double dh_dx = 
#if USE_RAYDIFF == 1
    raydifferentials_compute_rd_h(p, R, s, e);
#else
    raydifferentials_compute_rd_h(p, 0, s, e);
#endif
  if(dh_dx == 0.0) goto error;

  // nothing to be done here. do early out after compute_rd_h to init dh_dx
  if(e - s <= 1) return dh_dx;

  // call into camera to get precise per-pixel ray differential
  float rd_i[3], rd_j[3];
  if(raydifferentials_v1(p, 1.0, 1.0, rd_i, rd_j))
  { // mostly jumped out the viewing frustum
    stats->no_throughput++;
    return 0.0;
  }
  // express rd_i and rd_j not in world space but in tangent space of v1 (actually only 2d now):
  const float rduv_i[3] = {dotproduct(rd_i, p->v[1].diffgeo.dpdu), dotproduct(rd_i, p->v[1].diffgeo.dpdv), 0.0f};
  const float rduv_j[3] = {dotproduct(rd_j, p->v[1].diffgeo.dpdu), dotproduct(rd_j, p->v[1].diffgeo.dpdv), 0.0f};

  // transform rd_i,j from x[1] to h[i]
  const int camera_end = (s == 0) ? e : s;
  _halfvec_compute_stepsizes(p, R, rd_u, rd_v, rduv_i, rduv_j, camera_end);
  return dh_dx;
error:
  stats->raydiff_compute_failed ++;
  return 0.0;
}

// return dh_dx
static inline double halfvec_precache(
    halfvec_stats_t *stats,
    path_t *p,
    const int s,
    const int e,
    float *R,
    float *rd_u,
    float *rd_v)
{
  p->v[s].diffgeo.type = s_pinned_position;
  for(int i=s+1;i<=e;i++)
    p->v[i].diffgeo.type = s_free;
  // compute abc matrices:
  manifold_compute_tangents(p, s, e);
  return halfvec_compute_raydifferentials(stats, p, s, e, R, rd_u, rd_v);
}

static inline float halfvec_reverse_check(
    halfvec_stats_t *d,
    const path_t *curr,
    const path_t *tent,
    const int s,
    const int e)
{
  // check reversibility of markov chain, going from tent to curr.
  path_t reverse = *tent;
  float h[(e-s+1)*3];
  reverse.lambda = curr->lambda;
  reverse.time = curr->time;
  if(s == 0)
  { // mutated from camera, potentially changed point on lens:
    reverse.sensor.aperture_x = curr->sensor.aperture_x;
    reverse.sensor.aperture_y = curr->sensor.aperture_y;
    reverse.sensor.aperture_set = curr->sensor.aperture_set;
  }
  for(int v=s;v<e;v++)
  { // need to copy over iors, these depend on wavelength..
    reverse.e[v].vol.ior = curr->e[v].vol.ior;
    reverse.v[v].interior.ior = curr->v[v].interior.ior;
  }

  // both sub-path and aperture mutation need to match position of old starting point
  reverse.v[s] = curr->v[s];
  // recompute derivatives
  reverse.v[s].diffgeo.type = s_pinned_position;
  for(int i=s+1;i<=e;i++)
    reverse.v[i].diffgeo.type = s_free;
  if(manifold_compute_tangents(&reverse, s, e) == 0.0) goto error;

  // copy over target half vector constraints:
  h[0] = h[1] = h[2] = h[3*(e-s)] = h[3*(e-s)+1] = h[3*(e-s)+2] = 0.0f;
  for(int i=1;i<e-s;i++)
  {
    h[3*i+0] = curr->v[s+i].diffgeo.h[0];
    h[3*i+1] = curr->v[s+i].diffgeo.h[1];
    h[3*i+2] = curr->v[s+i].diffgeo.h[2];
  }

  const float pdf_reverse = halfvec_to_worldspace(d, &reverse, h, s, e);
  if(pdf_reverse == 0.0f) goto error;
  // also need to check world space positions:
  for(int i=s+1;i<e;i++)
  {
    float dist[3];
    float quant = 1.0f;
    for(int k=0;k<3;k++) quant = MAX(MAX(fabsf(tent->v[i].hit.x[k]), fabsf(curr->v[i].hit.x[k])), quant);
    for(int k=0;k<3;k++) dist[k] = reverse.v[i].hit.x[k] - curr->v[i].hit.x[k];
    const float len = sqrtf(dotproduct(dist, dist));
    if(len > quant * HALFVEC_REL_SPATIAL_EPS) goto error;
  }
  return pdf_reverse;
error:
  d->reverse_check_failed ++;
  return 0.0f;
}

// mutate distance of given edge e by (mirrored) gaussian with given sigma
static inline float halfvec_perturb_distance(const path_t *p, const int e, const float sigma)
{
  float g0, g1;
  const int tid = common_get_threadid();
  const float r0 = points_rand(rt.points, tid);
  const float r1 = points_rand(rt.points, tid);
  sample_gaussian(r0, r1, &g0, &g1);
  float dist = p->e[e].dist + sigma * g0;
  if(dist < 0)
    dist = -dist;
  dist = MAX(rt.epsilon, dist);
  return dist;
}

static inline float halfvec_pdf_perturb_distance(const path_t *p, const int e, const float sigma, const float dist)
{
  const double d0 =  p->e[e].dist - dist;
  const double d1 = -p->e[e].dist - dist;
  return 1.0f/(sigma*sqrtf(2.0f*M_PI)) *
    (exp(-d0*d0/(2.0f*sigma*sigma)) + exp(-d1*d1/(2.0f*sigma*sigma)));
}

static inline int _halfvec_perturb_internal(
    halfvec_stats_t *d,  // stats
    const path_t *curr,  // current path
    path_t *tent,        // tentative path, to be perturbed
    const int s,         // fixed vertex at the start of the glossy chain
    const int e,         // fixed vertex at the end of the glossy chain
    float *h,            // storage for the half vectors: 3*(e-s+1) floats
    const float *R,      // ray differential alignment matrix and main axes of curr path
    const float *rd_u,
    const float *rd_v)
{
  assert(e-s > 0);
  assert(s >= 0);
  assert(e < PATHSPACE_MAX_VERTS);

  const int tid = common_get_threadid();

  h[0] = h[1] = h[2] = h[3*(e-s)] = h[3*(e-s)+1] = h[3*(e-s)+2] = 0.0f;

  if(s == 0)
  {
    *tent = *curr; // wasteful, but whatever.

    // mutate wavelength a bit. this only needs a re-compute tangents if transmission happens somewhere on the path
    tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
    // TODO: try what happens if you mutate time, too (check with camera_mutate_aperture which camera time is used, tent or curr!)
    // TODO: need to just transport all vertices along their time derivatives (we have mostly linear motion)
    // TODO: and let mf fix it again.

    // mutate point on aperture, after mutating outgoing direction (pixel)
    // this works kind of well because we're offsetting the first point
    // only a little bit, the half vector to worldspace iteration will fix it.
    // XXX also do it on the light source!
    // XXX even if recomputing tangents, this causes a drop from 70-80 -- 50-60% acceptance!
    view_cam_mutate_aperture(tent, points_rand(rt.points, tid), points_rand(rt.points, tid), 0.2f);
    if(tent->v[0].mode == s_absorb) return 1; // sampling failed (stepped off the aperture etc)

    // only necessary if we mutate time, wavelength, or point on the camera or light.
    // we assume the caller takes care of his own business and hands us correct derivatives for tent.
    // so we only repair them in case we corrupted them by changing the path ourselves:
    if(manifold_compute_tangents(tent, s, e) <= 0.0f)
    {
      tent->length = 0;
      return 1;
    }
  }
  // else we assume the path was already inited and the end points were already inited or perturbed from the outside,
  // also constraint derivatives and transfer matrices need to be inited

  // mutate halfvectors
  for(int i=s+1;i<e;i++)
  {
    // sample gaussians based on roughness and ray diffs
    if(curr->v[i].mode & s_specular)
    {
      // force half vector to some precise (001) to avoid drift
      h[3*(i-s)+0] = 0.0f;
      h[3*(i-s)+1] = 0.0f;
      continue;
    }

    const float *currh = curr->v[i].diffgeo.h;
    float *tenth = h + 3*(i-s);
    float g0, g1, g2, g3, dh[3];
    sample_gaussian(points_rand(rt.points, tid), points_rand(rt.points, tid), &g0, &g1);
    sample_gaussian(points_rand(rt.points, tid), points_rand(rt.points, tid), &g2, &g3);

    const float s_bsdf = _halfvec_bsdf_stepsize(curr, i);
    const float s_dist = _halfvec_dist_stepsize(curr, i);
    const float s[3] = {s_bsdf, s_bsdf, s_dist};
    g0 *= rd_u[i];
    g1 *= rd_v[i];
    const float gv[3] = {g0, g1, g2};
    // tenth = currh + S * R * g
    mat3_mulv(R + 9*i, gv, dh);
    for(int k=0;k<3;k++) dh[k] *= s[k];
    for(int k=0;k<3;k++) tenth[k] = currh[k] + dh[k];
    if(!(curr->v[i].mode & s_volume)) tenth[2] = 0.0f;
    else if(tenth[2] < 0.0) tenth[2] = -tenth[2];
  }
  return 0;
}

// perturb the half vectors just the same as for HSLT, but instead
// of the newtonian solve, do a single predict/project step.
// return != 0 if any error ocurred
static inline double _halfvec_pdf_perturb_internal( halfvec_stats_t *stats, const path_t *curr, const path_t *tent, const int s, const int e, const float *h,
    const float *R, const float *rd_u, const float *rd_v);

static inline int halfvec_perturb_single(
    halfvec_stats_t *stats,  // stats
    const path_t *curr,      // current path
    path_t *tent,            // tentative path, to be perturbed
    const int s,             // fixed vertex at the start of the glossy chain
    const int e,             // fixed vertex at the end of the glossy chain
    const float *R,
    const float *rd_u,
    const float *rd_v)
{
  float h[3*(e-s+1)];
  if(_halfvec_perturb_internal(stats, curr, tent, s, e, h, R, rd_u, rd_v)) return 1;
  for(int i=1;i<e-s;i++)
  {
    // scale down by a bit, pretend we're proper HMC doing derivatives of path measurement contribution (bias towards h=0)
    // const float s_bsdf = _halfvec_bsdf_stepsize(curr, s+i);
    // const float hl = sqrtf(h[3*i+0]*h[3*i+0] + h[3*i+1]*h[3*i+1]);
    // const float scale = hl > s_bsdf ? (hl - s_bsdf/2.0)/hl : 1.0f;
    // const float scale = 0.8;
    // h[3*i+0] *= scale;
    // h[3*i+1] *= scale;
    h[3*i+0] -= curr->v[s+i].diffgeo.h[0];
    h[3*i+1] -= curr->v[s+i].diffgeo.h[1];
    h[3*i+2] -= curr->v[s+i].diffgeo.h[2];
  }
  if(manifold_map_h_to_x(curr, tent, h, 1.0, s, e))
  {
    tent->length = 0;
    return 1;
  }

  for(int v=s+1;v<e;v++)
  { // more stable projection scheme: first along normal and then check visibility via path_project()
    // DEBUG is this ortho normal?
    assert(fabsf(dotproduct(curr->v[v].diffgeo.dpdu, curr->v[v].diffgeo.dpdu) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(curr->v[v].diffgeo.dpdv, curr->v[v].diffgeo.dpdv) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(curr->v[v].hit.gn, curr->v[v].hit.gn) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(curr->v[v].diffgeo.dpdu, curr->v[v].diffgeo.dpdv)) < 1e-4f);
    assert(fabsf(dotproduct(curr->v[v].diffgeo.dpdu, curr->v[v].hit.gn)) < 1e-4f);
    assert(fabsf(dotproduct(curr->v[v].diffgeo.dpdv, curr->v[v].hit.gn)) < 1e-4f);

    // volume vertices aren't projected:
    if(curr->v[v].mode & s_volume) continue;

    ray_t ray;
    hit_t hit;
    ray.time = tent->time;
    ray.ignore = INVALID_PRIMID;
    for(int i=0;i<3;i++) ray.pos[i] = tent->v[v].hit.x[i];
    // project along normal direction, only need to project the point directly
    // to dpdu/dpdv for reverse walk, no distance computation for better accuracy.
    for(int i=0;i<3;i++) ray.dir[i] = curr->v[v].hit.gn[i];
    // TODO: could save a bit of compute by getting duv[] from map_h_to_x above.
    float dx[3], duv[2];
    for(int i=0;i<3;i++) dx[i] = tent->v[v].hit.x[i] - curr->v[v].hit.x[i];
    duv[0] = dotproduct(dx, curr->v[v].diffgeo.dpdu);
    duv[1] = dotproduct(dx, curr->v[v].diffgeo.dpdv);
    float max_dist = MAX(1e-4f, 4.0f*sqrtf(duv[0]*duv[0] + duv[1]*duv[1]));
    hit.prim = INVALID_PRIMID;
    hit.dist = max_dist;
    ray.min_dist = -max_dist;
    accel_closest(rt.accel, &ray, &hit, 0.0f);
    // path space things will be inited below in path_project (also checks visibility)
    if(hit.dist < max_dist && !primid_invalid(hit.prim))
      for(int i=0;i<3;i++) tent->v[v].hit.x[i] = ray.pos[i] + hit.dist * ray.dir[i];
    else
    {
      tent->length = v;
      stats->project_failed++;
      return 1;
    }
  }

  // fix up path vertex / edge data
  for(int v=s+1;v<=e;v++)
  { // project, trace towards position at v[v]
    float x[3] = {0.0}; // avoid jitter and set directly to position we projected above.
    if(v != e) for(int i=0;i<3;i++) x[i] = tent->v[v].hit.x[i];
    // this will also update the sensor struct/pixel coordinate if needed.
    tent->v[v].mode = curr->v[v].mode;
    tent->e[v].transmittance = 0.0f;
    if(path_project(tent, v, s_propagate_mutate) ||
        (tent->v[v-1].mode != curr->v[v-1].mode) ||
        (tent->v[v].flags  != curr->v[v].flags) ||
        (tent->v[v].hit.shader != curr->v[v].hit.shader) ||
        (tent->v[v].interior.shader != curr->v[v].interior.shader) ||
        (primid_invalid(tent->v[v].hit.prim) != primid_invalid(curr->v[v].hit.prim)))
    {
      tent->length = v;
      stats->project_failed++;
      return 1;
    }
    if(v != e)
    {
      const float eps = 1e-4f*fmaxf(fmaxf(.5f, fabsf(x[0])), fmaxf(fabsf(x[1]), fabsf(x[2])));
      for(int i=0;i<3;i++) if(fabsf(tent->v[v].hit.x[i] - x[i]) > eps)
      {
        tent->length = v;
        stats->project_failed++;
        return 1;
      }
      else tent->v[v].hit.x[i] = x[i];
    }
  }
  return 0;
}

// perturbs the tentative path's half vectors and transforms that back to
// worldspace.  if you're not perturbing from s=0, you should still initialize
// the ray differential matrices in R[s..e] and eigenvalues
// rd_uv[s..e], to guide anisotropic mutation.
static inline float halfvec_perturb(
    halfvec_stats_t *d,  // stats
    const path_t *curr,  // current path
    path_t *tent,        // tentative path, to be perturbed
    const int s,         // fixed vertex at the start of the glossy chain
    const int e,         // fixed vertex at the end of the glossy chain
    const float *R,      // ray differential matrices
    const float *rd_u,   // eigenvalues of curr path
    const float *rd_v)
{
  float h[3*(e-s+1)];
  if(_halfvec_perturb_internal(d, curr, tent, s, e, h, R, rd_u, rd_v)) return 0.0f;
  // convert half vectors to world space
  const float pdf = halfvec_to_worldspace(d, tent, h, s, e);
  if(pdf == 0.0f) return 0.0f; // detailed stats recorded downstream

  return pdf;
}

// returns a half vector space pdf
static inline double _halfvec_pdf_perturb_internal(
    halfvec_stats_t *stats,
    const path_t *curr,
    const path_t *tent,
    const int s,
    const int e,
    const float *h,
    const float *R,     // ray differential matric + eigenvalues of curr path
    const float *rd_u,
    const float *rd_v)
{
  double pdf = 1.0f;

  if(s == 0)
  {
    pdf *= view_cam_pdf_mutate_aperture(curr, tent, 0.2f);
    pdf *= spectrum_pdf_mutate(curr->lambda, tent->lambda);
    assert(pdf > 0.0);
  }

  for(int i=s+1;i<e;i++)
  {
    if(curr->v[i].mode & s_specular) continue;
    const float *tenth = (h ? h+3*(i-s) : tent->v[i].diffgeo.h), *currh = curr->v[i].diffgeo.h;
    float dp[3] = {tenth[0] - currh[0], tenth[1] - currh[1], tenth[2] - currh[2]};

    const float s_bsdf = _halfvec_bsdf_stepsize(curr, i);
    const float s_dist = _halfvec_dist_stepsize(curr, i);
    const float s[3] = {s_bsdf, s_bsdf, s_dist};
    float Rinv[9], d[3];
    mat3_transpose(R+9*i, Rinv);
    for(int k=0;k<3;k++) dp[k] /= s[k];
    mat3_mulv(Rinv, dp, d);
    const float v0 = rd_u[i];
    const float v1 = rd_v[i];
    const float v2 = 1.0f;
    if(curr->v[i].mode & s_volume)
    {
      // distance may have been mirrored, which is the sum of the gaussian
      // you would expect and the one around -currh[2] transformed through
      // the ray diff matrix:
      float dp2[3] = {tenth[0] - currh[0], tenth[1] - currh[1], tenth[2] + currh[2]}, d2[3];
      for(int k=0;k<3;k++) dp2[k] /= s[k];
      mat3_mulv(Rinv, dp2, d2);
      pdf *=   1.0f/(2.0f*M_PI*sqrtf(2.0f*M_PI) * v0*v1*v2) * exp(-0.5 * (d [0]*d [0]/(v0*v0) + d [1]*d [1]/(v1*v1) + d [2]*d [2]/(v2*v2)))
             + 1.0f/(2.0f*M_PI*sqrtf(2.0f*M_PI) * v0*v1*v2) * exp(-0.5 * (d2[0]*d2[0]/(v0*v0) + d2[1]*d2[1]/(v1*v1) + d2[2]*d2[2]/(v2*v2)));
    }
    else
      pdf *= 1.0f/(2.0f*M_PI*v0*v1) * exp(-0.5 * (d[0]*d[0]/(v0*v0) + d[1]*d[1]/(v1*v1)));
    // jacobian of scale (cancels out for all situations but textured roughness)
    pdf /= s[0]*s[1]*s[2];
  }
  return pdf;
}

static inline double halfvec_pdf_perturb(
    halfvec_stats_t *stats,
    const path_t *curr,
    const path_t *tent,
    const int s,
    const int e,
    const float *curr_R,
    const float *curr_rd_u,
    const float *curr_rd_v)
{
  return _halfvec_pdf_perturb_internal(stats, curr, tent, s, e, 0, curr_R, curr_rd_u, curr_rd_v);
}

// returns a half vector space pdf
static inline double halfvec_pdf_perturb_single(
    halfvec_stats_t *stats,
    const path_t *curr,
    const path_t *tent,
    const int s,
    const int e,
    const int reverse,
    const float *curr_R,
    const float *curr_rd_u,
    const float *curr_rd_v)
{
  double pdf_ratio = 1.0;

  for(int v=s+1;v<=e;v++)
  { // sanity checks corresponding to projection:
    if((tent->v[v-1].mode != curr->v[v-1].mode) ||
       (tent->v[v].flags  != curr->v[v].flags) ||
       (primid_invalid(tent->v[v].hit.prim) != primid_invalid(curr->v[v].hit.prim)))
      return 0.0;
  }
  // pdf of unprojection needs to account for ratio of geometry terms:
  float dx[3*PATHSPACE_MAX_VERTS] = {0.0};
  for(int v=s+1;v<e;v++)
  { // unproject tent->v[v] to tangent space of curr->v[v]
    float x[3] = {0.0f};
    for(int k=0;k<3;k++) x[k] = tent->v[v].hit.x[k] - curr->v[v].hit.x[k];
    // normal projection leaves dpdu/dpdv dotproducts unchanged:
    dx[3*(v-s)+0] = dotproduct(curr->v[v].diffgeo.dpdu, x);
    dx[3*(v-s)+1] = dotproduct(curr->v[v].diffgeo.dpdv, x);
    dx[3*(v-s)+2] = 0.0f;
    const float distn = dotproduct(curr->v[v].hit.gn, x);

    if(curr->v[v].mode & s_volume)
    { // volumes don't project, but have a 3d offset:
      dx[3*(v-s)+2] = distn;
      continue;
    }

    // ortho projection along normal incurs the following cosine correction factor:
    pdf_ratio *= path_lambert(tent, v, curr->v[v].hit.gn);

    if(reverse)
    { // shoot ray to make sure this was actually the closest point, but only for reverse pdf:
      ray_t ray;
      hit_t hit;
      ray.time = tent->time;
      ray.ignore = INVALID_PRIMID;
      for(int k=0;k<3;k++) ray.pos[k] = curr->v[v].hit.x[k] + dx[3*(v-s)+0] * curr->v[v].diffgeo.dpdu[k] + dx[3*(v-s)+1] * curr->v[v].diffgeo.dpdv[k];
      // for(int k=0;k<3;k++) ray.pos[k] = tent->v[v].hit.x[k] - distn * curr->v[v].hit.gn[k];
      for(int i=0;i<3;i++) ray.dir[i] = curr->v[v].hit.gn[i];
      float max_dist = MAX(1e-4f, 4.0f*sqrtf(dx[3*(v-s)+0]*dx[3*(v-s)+0] + dx[3*(v-s)+1]*dx[3*(v-s)+1]));
      if(fabsf(distn) > max_dist + 1e-4f) return 0.0; // too far, couldn't be our point.
      hit.prim = INVALID_PRIMID;
      hit.dist = max_dist;
      ray.min_dist = -max_dist;
      accel_closest(rt.accel, &ray, &hit, 0.0f);
      if(hit.dist >= max_dist || primid_invalid(hit.prim))
        return 0.0;
      else if(fabsf(distn - hit.dist) > 1e-4f) return 0.0;
    }
  }
  float h[3*PATHSPACE_MAX_VERTS];
  for(int v=s+1;v<e;v++)
  { // dh = map x to h (apply reverse matrix, i.e. without inversion)
    float a[3] = {0}, b[3], c[3] = {0};
    mat3_mulv(curr->v[v].diffgeo.a, dx+3*(v-s)-3, a);
    mat3_mulv(curr->v[v].diffgeo.b, dx+3*(v-s),   b);
    mat3_mulv(curr->v[v].diffgeo.c, dx+3*(v-s)+3, c);
    for(int k=0;k<3;k++)
      h[3*(v-s)+k] = curr->v[v].diffgeo.h[k] + a[k] + b[k] + c[k];
    // multiply up by a bit, pretend we're real HMC (see above)
    // const float s_bsdf = _halfvec_bsdf_stepsize(curr, v);
    // const float hl = sqrtf(h[3*(v-s)+0]*h[3*(v-s)+0] + h[3*(v-s)+1]*h[3*(v-s)+1]);
    // // this uses hl * x to compute scale, i.e. the unscaled variant
    // // XXX i think i got these upside down!
    // // XXX also the hl > ? test is now broken!
    // const float scale = hl > s_bsdf ? (hl + s_bsdf/2.0)/hl : 1.0f;
    // const float scale = 0.8;
    // h[3*(v-s)+0] /= scale;
    // h[3*(v-s)+1] /= scale;
    // pdf_ratio *= scale * scale;
  }

  double pdf_h = _halfvec_pdf_perturb_internal(stats, curr, tent, s, e, h, curr_R, curr_rd_u, curr_rd_v);

  // to convert to vertex area measure (or use halfvec measurement)
  // you need to multiply curr dh/dx, not tent dh/dx as would be the case
  // for the full half vector iteration, converging to h_t=h_p.
  return pdf_ratio * pdf_h;
}

// reflect the incoming direction at the give vertex (p->e[v].omega) around the
// half vector there, and initialize the outgoing direction at the next segment
// (p->e[v+1].omega). might fail if total internal reflection occurs.
static inline int halfvec_reflect(path_t *p, int v)
{
  if((p->v[v].mode & s_volume) || (p->v[v].mode & s_fiber))
  { // project to riemann sphere:
    float h[3] = {p->v[v].diffgeo.h[0], p->v[v].diffgeo.h[1], 1.0};
    const float z = 2.0f/(h[0]*h[0] + h[1]*h[1] + 1.0f);
    h[0] *= z;
    h[1] *= z;
    h[2] = z - 1.0f;
    for(int k=0;k<3;k++)
      p->e[v+1].omega[k] =
        p->v[v].diffgeo.dpdu[k] * h[0] +
        p->v[v].diffgeo.dpdv[k] * h[1] +
        p->v[v].hit.gn[k] * h[2];
    return 0;
  }

  float h[3], H[3]; // reconstruct world space half vector
  // plane/plane half vector lives in ortho normal basis hit.{a,b,n}
  // (and usually points into hemisphere of optically thinner medium, but we don't care)
  h[0] = p->v[v].diffgeo.h[0];
  h[1] = p->v[v].diffgeo.h[1];
  h[2] = 1.0f/sqrtf(1.0 + h[0]*h[0] + h[1]*h[1]);
  h[0] *= h[2];
  h[1] *= h[2];
  // h points into random hemisphere (always aligns with n), but sign of H
  // doesn't matter further down in the computation
  for(int k=0;k<3;k++)
    H[k] = p->v[v].hit.a[k] * h[0] +
      p->v[v].hit.b[k] * h[1] +
      p->v[v].hit.n[k] * h[2];

  // now update outgoing direction accordingly
  if(p->v[v].mode & s_reflect)
  {
    const float dot = dotproduct(p->e[v].omega, H);
    for(int k=0;k<3;k++)
      p->e[v+1].omega[k] = p->e[v].omega[k] - 2.0f*H[k] * dot;
  }
  else if(p->v[v].mode & s_transmit)
  {
    const float eta_ratio = mf(path_eta_ratio(p, v), 0); // = n1/n2;
    if(fabsf(eta_ratio - 1.0f) > 1e-5f)
    {
      float cosr = dotproduct(H, p->e[v].omega);
      if(cosr > 0.0f) // need to flip H to point into R hemisphere
        for(int k=0;k<3;k++) H[k] = -H[k];
      else cosr = -cosr;
      if(eta_ratio < 0.0f)
        return 1; // volume nesting broken.
      const float cost2 = 1.0f - eta_ratio*eta_ratio * (1.0f - cosr*cosr);
      if(cost2 <= 0.0f)
        return 1; // total internal reflection
      const float cost = sqrtf(cost2);
      const float f = eta_ratio*cosr - cost;
      for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k]*eta_ratio + f * H[k];
      normalise(p->e[v+1].omega);
    }
    else
    {
      for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k];
    }
  }
  return 0;
}


#undef HALFVEC_BSDF_STEP
#undef HALFVEC_BECKMANN_MIN
#undef HALFVEC_BECKMANN_MAX

#undef USE_RAYDIFF

#endif

