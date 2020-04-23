#pragma once

#include "points.h"
#include "pathspace/halfvec.h"
#include "pathspace/multichain.h"

#include "shader.h"
#include "camera.h"

static inline void hslt_get_b_cdf(path_t *p, int min_b, float *cdf)
{
  // actually sampling hair as a breakup point is fine,
  // just the chains shouldn't span hair.
  // // looking at hair?
  // int first_valid_b = min_b;
  int has_hair = 0;
  // manifold walks do so far not work with hair scattering :(
  for(int k=1;k<p->length;k++) if(p->v[k].mode & s_fiber)
  {
    has_hair = 1;
    break;
  }
      // first_valid_b = k; // b on a hair vertex is fine, it's halfvector will not be perturbed explicitly

  // by roughness only
  float sum_rough = 0.0f;
  if(has_hair) cdf[0] = 0.0f;
  else cdf[0] = 0.1f;      // full halfvec
  cdf[p->length-1] = 0.1f; // mull multi chain
  for(int k=1;k<p->length-1;k++)
  {
    int indexmatched = ((p->v[k].mode & s_transmit) && fabsf(path_eta_ratio(p, k) - 1.0f) < 1e-3f);
    float step = p->v[k].shading.roughness;
    const float g = p->v[k].interior.mean_cos;
    if(p->v[k].mode & s_volume) step = sqrtf(1.0f - g*g);
    sum_rough += (cdf[k] = (indexmatched ? 0 : step));
    // prefer hair points that can be connected to non-hairs:
    if((p->v[k].mode & s_fiber) &&
      !(p->v[k+1].mode & s_fiber)) cdf[k] *= 10.0f;
    //cdf[k] *= powf(4.0f, -k); // prefer early connections for inreased speed
  }

  if(sum_rough == 0.0f) cdf[0] = 0.0f; // no pure halfvec possible.

  // make cdf:
  for(int k=1;k<p->length;k++)
    cdf[k] += cdf[k-1];
  for(int k=0;k<p->length;k++)
    cdf[k] /= cdf[p->length-1];
}

static inline void hslt_get_c_cdf(const path_t *p, int b, float *cdf)
{
  // don't want to place c at vert 1, ever.
  b = MAX(1, b);
  for(int k=0;k<=b;k++) cdf[k] = 0.0f;
  // get c ~ 1/G(c-1, c)
  for(int k=b+1;k<p->length-1;k++)
  {
    if(p->v[k].mode & s_fiber)
    { // don't want a half vec chain to span across fiber scattering
      // events, we don't have good derivatives for it (nor do i
      // expect visibility of such connections to play nicely)
      cdf[k++] = 1e10f;
      for(;k<p->length;k++) cdf[k] = 0.0f;
      break;
    }
    // TODO: encourage c such that crazy diffgeo between c and b doesn't happen!
    int indexmatched = ((p->v[k].mode & s_transmit) && fabsf(path_eta_ratio(p, k) - 1.0f) < 1e-3f);
    float step = p->v[k].shading.roughness;
    const float g = p->v[k].interior.mean_cos;
    if(p->v[k].mode & s_volume) step = sqrtf(1.0f - g*g);
    cdf[k] = p->e[k].dist * p->e[k].dist * (indexmatched ? 0 : step);
  }
  cdf[p->length-1] = p->e[p->length-1].dist * p->e[p->length-1].dist;

  for(int k=b+1;k<p->length;k++)
    cdf[k] *= powf(2.0f, -(k-b-1)); // prefer early connections for inreased speed

  // make cdf:
  for(int k=b+2;k<p->length;k++)
    cdf[k] += cdf[k-1];
  for(int k=b+1;k<p->length-1;k++)
    cdf[k] /= cdf[p->length-1];
  cdf[p->length-1] = 1.0f;

  // XXX DEBUG force c = b+1
  // memset(cdf, 0, sizeof(float)*p->length);
  // for(int k=b+1;k<p->length+1;k++)
  //   cdf[k] = 1.0f;
}

static inline float hslt_perturb(
    path_t *curr,       // constant current sample, needs inited half vectors
    path_t *tent,       // tentative sample with tent->sensor inited to the desired sensor offset
    const int b,        // breakup vertex
    const int c,        // last vertex, where connecting back to base path
    const int single,   // do single step instead of newtonian walk
    halfvec_stats_t *stats)
{
  // copy sensor struct
  sensor_t sensor = tent->sensor;

  // first very wastefully copy the whole path
  // TODO: really not necessary, could only do vertex by vertex (multichain inits them anyhow)
  // TODO: maybe do this, but then construct multichain on different, empty path
  // and use this to converge the halfvector space part in the same index range?
  // XXX actually the vmlt_hslt.h already copied this and then updated the sensor struct!
  *tent = *curr;
  tent->sensor = sensor; // restore
  tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, common_get_threadid()), 0);
  tent->time = curr->time; // <= TODO? or use time perturbation instead.

  assert(!(tent->v[b].mode & s_specular));
  assert(!(tent->v[c].mode & s_specular));

  // track pdf ratio T(t->c)/T(c->t) in vertex area measure (used for camera vertex)
  double T = multichain_perturb(curr, tent, 0, b);
  if(T <= 0.0) return 0.0f;

  // done doing multichain perturbation.
  // compute offset from current to tentative index counting:
  // tent->v[v+to] corresponds to curr->v[v] for v >= b.
  const int to = tent->length - curr->length;

  // update segment (b..c) by half vector perturbation
  double pdf_fwd = 1.0, pdf_rev = 1.0;
  path_t hvec_storage;
  path_t *hvec = tent;
  if(c-b > 1)
  { // work to do, halfvector mutation
    if(to)
    {
      // potentially port half vectors to different index scheme.  fastest path
      // for a quick hack is another copy of a path.. :( we only do this in
      // case it's really required (the to offset is > 0)
      hvec = &hvec_storage;
      *hvec = *curr;
      hvec->time = tent->time;
      hvec->lambda = tent->lambda;
      // copy over mutated point b and incoming segment (medium/ior) from
      // multichain perturbation
      hvec->v[b] = tent->v[b+to];
      hvec->e[b] = tent->e[b+to];
    }

    if(single)
    {
      float curr_R[9*PATHSPACE_MAX_VERTS];
      float curr_rd_u[PATHSPACE_MAX_VERTS];
      float curr_rd_v[PATHSPACE_MAX_VERTS];
      double curr_dh_dx = halfvec_precache(stats, curr, b, c, curr_R, curr_rd_u, curr_rd_v);
      if(curr_dh_dx == 0.0) return 0.0f;

      // create a new path by perturbing half vectors:
      if(halfvec_perturb_single(stats, curr, hvec, b, c, curr_R, curr_rd_u, curr_rd_v))
      {
        stats->perturb_failed++;
        return 0.0f;
      }

      // TODO: smooth vertices
      float hvec_R[9*PATHSPACE_MAX_VERTS];
      float hvec_rd_u[PATHSPACE_MAX_VERTS];
      float hvec_rd_v[PATHSPACE_MAX_VERTS];
      double hvec_dh_dx = halfvec_precache(stats, hvec, b, c, hvec_R, hvec_rd_u, hvec_rd_v);
      if(hvec_dh_dx == 0.0) return 0.0f;

      // compute transition probabilities:
      pdf_fwd = curr_dh_dx * halfvec_pdf_perturb_single(stats, curr, hvec, b, c, 0, curr_R, curr_rd_u, curr_rd_v);
      if(!(pdf_fwd > 0.0f))
      {
        stats->no_throughput++;
        return 0.0f;
      }

      // computing the reverse pdf requires ray differentials to be initialized on tent
      pdf_rev = hvec_dh_dx * halfvec_pdf_perturb_single(stats, hvec, curr, b, c, 1, hvec_R, hvec_rd_u, hvec_rd_v);
      if(!(pdf_rev > 0.0f))
      {
        stats->reverse_check_failed++;
        return 0.0f;
      }
    }
    else
    {
      // re-init constraint derivatives for sub-path (b..c)
      // also doing this for current sample already.
      // we need the derivatives for the transfer matrix below,
      // but this also inits the half vectors:
      // TODO: could be optimized by just slapping over the a b c matrices for v[b]!
      // TODO: the rest should already be inited by the full path tangent compute above..
      curr->v[b].diffgeo.type = hvec->v[b].diffgeo.type = s_pinned_position;
      for(int i=b+1;i<=c;i++)
        curr->v[i].diffgeo.type = hvec->v[i].diffgeo.type = s_free;
      if(manifold_compute_tangents(hvec, b, c) == 0.0 ||
          manifold_compute_tangents(curr, b, c) == 0.0)
        return 0.0f;

      // create a new path by perturbing half vectors:
      float curr_R[9*PATHSPACE_MAX_VERTS];
      float curr_rd_u[PATHSPACE_MAX_VERTS];
      float curr_rd_v[PATHSPACE_MAX_VERTS];
      const double curr_dh_dx = halfvec_compute_raydifferentials(stats, curr, b, c, curr_R, curr_rd_u, curr_rd_v);
      if(curr_dh_dx == 0.0) return 0.0f;
      const float vol_pdf = halfvec_perturb(stats, curr, hvec, b, c, curr_R, curr_rd_u, curr_rd_v);
      if(vol_pdf <= 0.0f) return 0.0f;

      // check if we can actually walk back, to satisfy detailed balance
      const float vol_rpdf = halfvec_reverse_check(stats, curr, hvec, b, c);
      if(vol_rpdf <= 0.0f) return 0.0f;

      // compute transition probabilities:
      pdf_fwd = halfvec_pdf_perturb(stats, curr, hvec, b, c, curr_R, curr_rd_u, curr_rd_v);
      if(pdf_fwd <= 0.0f) return 0.0f;

      // compute transfer matrix for tentative sample:
      float hvec_R[9*PATHSPACE_MAX_VERTS];
      float hvec_rd_u[PATHSPACE_MAX_VERTS];
      float hvec_rd_v[PATHSPACE_MAX_VERTS];
      const double hvec_dh_dx = halfvec_compute_raydifferentials(stats, hvec, b, c, hvec_R, hvec_rd_u, hvec_rd_v);
      if(hvec_dh_dx == 0.0) return 0.0f;
      pdf_rev = halfvec_pdf_perturb(stats, hvec, curr, b, c, hvec_R, hvec_rd_u, hvec_rd_v);
      if(pdf_rev <= 0.0f) return 0.0f;
    }
    // construct final tentative path by copying half vector perturbed postfix (b+1..c) over
    if(to)
    {
      memcpy(tent->v + b+1+to, hvec->v + b+1, sizeof(vertex_t)*(c - b));
      memcpy(tent->e + b+1+to, hvec->e + b+1, sizeof(edge_t)*(c - b));
    }
  }
  else if(c > b)
  { // simple connection
    // quite terrible hack: don't test this for fiber/fiber connections
    // if(!((curr->v[b].mode & s_fiber) && (curr->v[c].mode & s_fiber)))
    {
      if(path_project(tent, c+to, s_propagate_mutate) ||
          (tent->v[c+to].flags != curr->v[c].flags) ||
          (tent->v[c+to].hit.shader != curr->v[c].hit.shader) ||
          (tent->v[c+to].interior.shader != curr->v[c].interior.shader) ||
          (primid_invalid(tent->v[c+to].hit.prim) != primid_invalid(curr->v[c].hit.prim))) 
        return 0.0f;
      // check whether we actually arrived at vertex c
      for(int k=0;k<3;k++)
        if(fabsf(tent->v[c+to].hit.x[k] - curr->v[c].hit.x[k]) > HALFVEC_REL_SPATIAL_EPS *
            MAX(MAX(fabsf(tent->v[c+to].hit.x[k]), fabsf(curr->v[c].hit.x[k])), 1.0))
          return 0.0f;
    }
#if 0
    else
    {
      for(int k=0;k<3;k++) tent->e[c+to].omega[k] = tent->v[c+to].hit.x[k] - tent->v[c+to-1].hit.x[k];
      tent->e[c+to].dist = sqrtf(dotproduct(tent->e[c+to].omega, tent->e[c+to].omega));
      for(int k=0;k<3;k++) tent->e[c+to].omega[k] /= tent->e[c+to].dist;
      // if(tent->e[c+to].dist < 0.001f) return 0.0f;
      tent->e[c+to].dist = MAX(tent->e[c+to].dist, 0.1f);
    }
#endif
  }
  // else c == b and no connection necessary

  // catch stupid special case that we bounced off a light source and the
  // measurement contribution actually comes from this sub-path and not from
  // the last vertex.
  if(lights_eval_vertex(tent, tent->length-1) <= 0.0) return 0.0f;

  // compute acceptance in vertex area measure:

  // need to evaluate measurement contributions from start to end of path,
  // as shader_brdf() will init the reflect/transmit flags which may
  // affect volume stack madness in path_ior_ratio().
  double f_tent = 1.0;
  double f_curr = 1.0;

  // eval connection point (bsdf at b) to init inner vertex flags
  if(b > 0 && b != c) f_tent *= shader_brdf(tent, b+to);
  if(b > 0 && b != c) f_curr *= shader_brdf(curr, b);
  if(single)
  {
    f_tent *= path_measurement_contribution_dx(hvec, b, c);
    f_curr *= path_measurement_contribution_dx(curr, b, c);
  }
  else if(c > b)
  { // degenerates to dx measurement for c = b+1:
    f_tent *= halfvec_measurement(hvec, b, c);
    f_curr *= halfvec_measurement(curr, b, c);
  }
  // same at c
  if(c != curr->length-1) f_tent *= shader_brdf(tent, c+to);
  if(c != curr->length-1) f_curr *= shader_brdf(curr, c);

  if(f_tent <= 0.0f || f_curr <= 0.0f)
    stats->no_throughput++;
  stats->mutations_proposed++;

  assert(f_tent <= 0.0 || path_measurement_contribution_dwp(tent, 0, tent->length-1) > 0);

  // compute acceptance as f->/T->  /  f<-/T<-
  return T * f_tent/pdf_fwd / f_curr*pdf_rev;
}
