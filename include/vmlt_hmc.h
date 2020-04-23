#pragma once

#include "pathspace/halfvec.h"
#include "tools/geo/texture.h"
#include "tools/geo/vdata.h"
#include "geo.h"
#include "geo/triangle.h"

typedef struct vmlt_hmc_t
{
  halfvec_stats_t *stats;
  vdata_t *normals, *footprint;
  texture_t *leadr;
}
vmlt_hmc_t;

void *hmc_init()
{
  vmlt_hmc_t *d = (vmlt_hmc_t *)malloc(sizeof(vmlt_hmc_t));
  memset(d, 0, sizeof(vmlt_hmc_t));
  d->stats = (halfvec_stats_t *)malloc(sizeof(halfvec_stats_t)*rt.num_threads);
  memset(d->stats, 0, rt.num_threads*sizeof(halfvec_stats_t));

  // get leadr data for the shapes:
  d->footprint = (vdata_t *)malloc(sizeof(vdata_t)*rt.prims->num_shapes);
  d->normals = (vdata_t *)malloc(sizeof(vdata_t)*rt.prims->num_shapes);
  d->leadr = (texture_t *)malloc(sizeof(texture_t)*rt.prims->num_shapes);
  for(int s=0;s<rt.prims->num_shapes;s++)
  {
    char filename[512];
    snprintf(filename, sizeof(filename), "%s/%s-leadr.tex", rt.searchpath, rt.prims->shape[s].name);
    if(!texture_load(d->leadr+s, filename, 1))
      fprintf(stderr, "[hmc] loaded leadr texture %s\n", filename);
    snprintf(filename, sizeof(filename), "%s/%s-leadr.vdata", rt.searchpath, rt.prims->shape[s].name);
    if(!vdata_map(d->normals+s, filename))
      fprintf(stderr, "[hmc] loaded base normals %s\n", filename);
    snprintf(filename, sizeof(filename), "%s/%s.vdata", rt.searchpath, rt.prims->shape[s].name);
    if(!vdata_map(d->footprint+s, filename))
      fprintf(stderr, "[hmc] loaded micropolygon uv footprints %s\n", filename);
    // stupid file format. should change to load this directly:
    d->normals[s].num_slots = 3;
    d->footprint[s].num_slots = 5;
  }

  return d;
}

static inline void hmc_stats_accum(halfvec_stats_t *accum, const vmlt_hmc_t *d)
{
  memset(accum, 0, sizeof(halfvec_stats_t));
  for(int k=0;k<rt.num_threads;k++)
  {
    accum->project_failed += d->stats[k].project_failed;
    accum->propagate_failed += d->stats[k].propagate_failed;
    accum->perturb_failed += d->stats[k].perturb_failed;
    accum->no_throughput += d->stats[k].no_throughput;
    accum->mutations += d->stats[k].mutations;
    accum->mutations_proposed += d->stats[k].mutations_proposed;
    accum->reverse_check_failed += d->stats[k].reverse_check_failed;
    accum->sum_iterations += d->stats[k].sum_iterations;
    accum->sum_successful_reduces += d->stats[k].sum_successful_reduces;
    accum->raydiff_compute_failed += d->stats[k].raydiff_compute_failed;
    accum->error_increased += d->stats[k].error_increased;
  }
}

void hmc_cleanup(void *data)
{
  vmlt_hmc_t *d = (vmlt_hmc_t *)data;
  for(int s=0;s<rt.prims->num_shapes;s++)
  {
    texture_cleanup(d->leadr+s);
    vdata_cleanup(d->normals+s);
    vdata_cleanup(d->footprint+s);
  }
  free(d->leadr);
  free(d->normals);
  free(d->stats);
  free(d);
}

float hmc_suitability(const path_t *p, void *data)
{
  // need a complete path with at least one bsdf vertex
  if(p->length <= 2) return 0.0f;
  // non specular only!
  for(int i=1;i<p->length-1;i++)
    if(((p->v[i].mode & s_volume) == 0) &&
        (p->v[i].shading.roughness == 0.0f)) return 0.0f;
  return 1.0f;
}

static inline void _hmc_get_ortho(
    const float *n,  // new normal
    float *dpdu,     // return dpdu but perpendicular to n
    float *dpdv,     // return dpdv but perp u and perp n.
    float *dndu,     // also orthogonalise normal derivative, if != 0
    float *dndv)
{
  // gram schmidt:
  const float ilu = 1.0f/sqrtf(dotproduct(dpdu, dpdu));
  for(int k=0;k<3;k++) dpdu[k] *= ilu;
  if(dndu) for(int k=0;k<3;k++) dndu[k] *= ilu;
  const float dotu = dotproduct(dpdu, dpdv);
  for(int k=0;k<3;k++) dpdv[k] -= dotu*dpdu[k];
  if(dndv) for(int k=0;k<3;k++) dndv[k] -= dotu*dndu[k];
  const float dotn = dotproduct(n, dpdv);
  for(int k=0;k<3;k++) dpdv[k] -= dotn*n[k];
  if(dndv) for(int k=0;k<3;k++) dndv[k] -= dotn*n[k];
  const float ilv = 1.0f/sqrtf(dotproduct(dpdv, dpdv));
  for(int k=0;k<3;k++) dpdv[k] *= ilv;
  if(dndv) for(int k=0;k<3;k++) dndv[k] *= ilv;
}

static inline void _hmc_lookup_vertex(
    const vmlt_hmc_t *d,  // global textures/vertex data
    const int vi,         // vertex index on primitive to lookup
    const path_t *p,      // path
    const int v,          // path vertex
    const int level,      // smoothness level (0 means don't do anything)
    float *v0, float *v1, float *v2,    // vertex positions in world space
    float *uv0, float *uv1, float *uv2, // corresponding uv coordinates in texture space
    float *lod_n,         // return vertex normal
    float *lod_rough)     // return roughness u,v,uv (in dpdu/dpdv space)
{
  const int s  = p->v[v].hit.prim.shapeid;
  const int mb = p->v[v].hit.prim.mb + 1;

  // get 5 leadr moments in uv space from texture.
  // do ewa texture lookup multiplied by 2^level
  const float scale = powf(2.0f, level);
  float fu, fv, d0[2], d1[2], leadr_E[5];
  const float *vd = vdata_get(d->footprint+s, vi);
  for(int k=0;k<2;k++) d0[k] = scale * vd[1+k];
  for(int k=0;k<2;k++) d1[k] = scale * vd[3+k];
  geo_decode_uv(((uint32_t*)vd)[0], &fu, &fv);
  texture_lookup_ewa(d->leadr+s, fu, fv, d0, d1, leadr_E);

  // get undisplaced base normal from vdata, for current time:
  float n[3], dpdu[3], dpdv[3];
  if(mb > 1)
  {
    const float *n0 = vdata_get(d->normals+s, mb*vi);
    const float *n1 = vdata_get(d->normals+s, mb*vi+1);
    for(int k=0;k<3;k++) n[k] = (1.0f-p->time)*n0[k] + p->time*n1[k];
    normalise(n);
  }
  else
  {
    const float *n0 = vdata_get(d->normals+s, vi);
    for(int k=0;k<3;k++) n[k] = n0[k];
  }
  // get dpdu dpdv around this normal
  for(int k=0;k<3;k++) dpdu[k] = p->v[v].diffgeo.dpdu[k];
  for(int k=0;k<3;k++) dpdv[k] = p->v[v].diffgeo.dpdv[k];
  _hmc_get_ortho(n, dpdu, dpdv, 0, 0);

  // set up transformation from uv to dpdu/dpdv/n:
  // for that, we need a shared space to transition between tangent frame
  // and uv coordinates. this will be a `barycentric' space that pretends v0 v1 v2 span
  // a triangle and our undisplaced vertex and all neighbours lie in this plane.
  float Mp[4], Muv[4]; // convert from barycentric to p and uv spaces
  float e0[3], e1[3];
  for(int k=0;k<3;k++) e0[k] = v2[k] - v0[k];
  for(int k=0;k<3;k++) e1[k] = v1[k] - v0[k];
  // rotate from `barycentric' to tangent frame dpdu/dpdv.
  // actually these are not in the same plane, but we pretend the
  // undisplaced triangle edges would be. the dot is the same.
  Mp[0] = dotproduct(dpdu, e0);
  Mp[1] = dotproduct(dpdu, e1);
  Mp[2] = dotproduct(dpdv, e0);
  Mp[3] = dotproduct(dpdv, e1);
  Muv[0] = uv2[0] - uv0[0];
  Muv[1] = uv1[0] - uv0[0];
  Muv[2] = uv2[1] - uv0[1];
  Muv[3] = uv1[1] - uv0[1];
  // uv = Muv * Mp^{-1} * tangent frame vector
  float M[4], iMp[4];
  mat2_invert(Mp, iMp);
  mat2_mul(Muv, iMp, M);
  // M = | du/dx du/dy |
  //     | dv/dx dv/dy |
  // convert leadr moments from uv to dpdu/dpdv space (see leadr paper eq (9) and (12)
  const float Ex = M[0] * leadr_E[0] + M[2] * leadr_E[1];
  const float Ey = M[1] * leadr_E[0] + M[3] * leadr_E[1];
  const float s2u = leadr_E[2] - leadr_E[0]*leadr_E[0];
  const float s2v = leadr_E[3] - leadr_E[1]*leadr_E[1];
  const float cuv = leadr_E[4] - leadr_E[0]*leadr_E[1];
  const float s2x = M[0]*M[0] * s2u + M[2]*M[2] * s2v + 2.0f*M[0]*M[2]*cuv;
  const float s2y = M[1]*M[1] * s2u + M[3]*M[3] * s2v + 2.0f*M[1]*M[3]*cuv;
  const float cxy = M[0]*M[1]*s2u + M[2]*M[3]*s2v + (M[0]*M[1]+M[2]*M[3])*cuv;
  // n = normalise (E(x), E(y), 1)
  lod_n[0] = Ex;
  lod_n[1] = Ey;
  lod_n[2] = 1.0f;
  normalise(lod_n);
  lod_rough[0] = s2x;
  lod_rough[1] = s2y;
  lod_rough[2] = cxy;
}

// overwrite vertex differential geometry and roughness:
static inline void _hmc_smooth_vertex(
    vmlt_hmc_t *d,
    path_t *p,        // path to alter
    const int v,      // path vertex to alter
    const int level)  // smoothness level (0 means don't do anything)
{
  if(level == 0) return; // now that was easy.

  const int s  = p->v[v].hit.prim.shapeid;
  const int vi = p->v[v].hit.prim.vi;

  if(!d->leadr[s].tex) return;
  if(!d->normals[s].data) return;
  if(!d->footprint[s].data) return;

  // do a lookup for all three adjacent vertices in the mesh:
  // vertex locations and uv coordinates
  float v0[3], v1[3], v2[3], uv0[2], uv1[2], uv2[2];
  // barycentric coordinates:
  float b0, b1;
  // vertex indices to do lookups on:
  int vi0 = rt.prims->shape[s].vtxidx[vi + 0].v, vi1;
  int vi2 = rt.prims->shape[s].vtxidx[vi + 2].v;

  geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 0, p->time, v0);
  geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 2, p->time, v2);
  geo_get_uv(rt.prims, p->v[v].hit.prim, 0, uv0);
  geo_get_uv(rt.prims, p->v[v].hit.prim, 2, uv2);

  if(p->v[v].hit.prim.vcnt == 3)
  {
    vi1 = rt.prims->shape[s].vtxidx[vi + 1].v,
    geo_get_uv(rt.prims, p->v[v].hit.prim, 1, uv1);
    geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 1, p->time, v1);
    b0 = p->v[v].hit.u;
    b1 = p->v[v].hit.v;
  }
  else if(p->v[v].hit.prim.vcnt == 4)
  {
    if(p->v[v].hit.v >= p->v[v].hit.u)
    {
      vi1 = rt.prims->shape[s].vtxidx[vi + 1].v,
      geo_get_uv(rt.prims, p->v[v].hit.prim, 1, uv1);
      geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 1, p->time, v1);
      b0 = p->v[v].hit.u;
      b1 = p->v[v].hit.v - p->v[v].hit.u;
    }
    else
    {
      vi1 = rt.prims->shape[s].vtxidx[vi + 3].v,
      geo_get_uv(rt.prims, p->v[v].hit.prim, 3, uv1);
      geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 3, p->time, v1);
      b0 = p->v[v].hit.u - p->v[v].hit.v;
      b1 = p->v[v].hit.v;
    }
  }
  else return; // no mesh, no smoothing.

  // get meso normals and roughnesses for three mesh vertices
  float n0[3], n1[3], n2[3], r0[3], r1[3], r2[3];
  _hmc_lookup_vertex(d, vi0, p, v, level, v0, v1, v2, uv0, uv1, uv2, n0, r0);
  _hmc_lookup_vertex(d, vi1, p, v, level, v0, v1, v2, uv0, uv1, uv2, n1, r1);
  _hmc_lookup_vertex(d, vi2, p, v, level, v0, v1, v2, uv0, uv1, uv2, n2, r2);

  // very importantly get dndu/dndv, first in barycentric triangle space:
  geo_tri_dnduv(0, 0, 0, n0, n1, n2, uv0, uv1, uv2, b0, b1, 0, 0, p->v[v].diffgeo.dndu, p->v[v].diffgeo.dndv);

  // interpolate hit->n from three surrounding normals
  const float b2 = 1.0f - b0 - b1;
  for(int k=0;k<3;k++) p->v[v].hit.n[k] = b0 * n2[k] + b1 * n1[k] + b2 * n0[k];
  normalise(p->v[v].hit.n);
  for(int k=0;k<3;k++) p->v[v].hit.gn[k] = p->v[v].hit.n[k];

  // update diffgeo.dpdu/dpdv once again to be ortho to hit.n
  _hmc_get_ortho(p->v[v].hit.n, p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv, p->v[v].diffgeo.dndu, p->v[v].diffgeo.dndv);
  
  // roughness:
  // TODO: do we really do this? roughness is computed around undisplaced normal
  //       in pre-pass! can we just rotate it to new tangent frame?
  // set up transformation from uv to dpdu/dpdv/n
  // Mp[0] = dotproduct(p->v[v].diffgeo.dpdu, e0);
  // Mp[1] = dotproduct(p->v[v].diffgeo.dpdu, e1);
  // Mp[2] = dotproduct(p->v[v].diffgeo.dpdv, e0);
  // Mp[3] = dotproduct(p->v[v].diffgeo.dpdv, e1);
  p->v[v].shading.roughness_type = s_rough_anisotropic_beckmann;
  // interpolate using barycentric coordinates
  // XXX this is wrong, should interpolate the E(x) moments istead.
  p->v[v].shading.roughness    = b0 * r2[0] + b1 * r1[0] + b2 * r0[0];
  p->v[v].shading.roughness_v  = b0 * r2[0] + b1 * r1[1] + b2 * r0[1];
  p->v[v].shading.roughness_uv = b0 * r2[0] + b1 * r1[1] + b2 * r0[1];
}


float hmc_mutate(
    path_t *curr,
    path_t *tent,
    void *data)
{
  vmlt_hmc_t *d = (vmlt_hmc_t *)data;
  halfvec_stats_t *stats = d->stats + common_get_threadid();
  // halfvec_stats_t dummy = {0};
  // if(rt.pointsampler->t[common_get_threadid()].num_rejects < 1000)
  // if(curr->length != 5)
    // stats = &dummy;
  stats->mutations++;

  // TODO: sample a smoothing level. one or one per vertex?

  // if curr wasn't created using a half vector mutation there's a good chance ray diffs will not be inited:
  // if(!curr->cache.inited) // actually we don't need the cache.f measurement contribution.
  // TODO: smooth vertices
  // scramble volume tangent frame orientation
  curr->tangent_frame_scrambling = 0.1f + points_rand(rt.points, common_get_threadid())*(0.9f-0.1f);
  for(int k=1;k<curr->length-1;k++) manifold_init(curr, k);
  float curr_R[9*PATHSPACE_MAX_VERTS];
  float curr_rd_u[PATHSPACE_MAX_VERTS];
  float curr_rd_v[PATHSPACE_MAX_VERTS];
  float tent_R[9*PATHSPACE_MAX_VERTS];
  float tent_rd_u[PATHSPACE_MAX_VERTS];
  float tent_rd_v[PATHSPACE_MAX_VERTS];
  double curr_dh_dx = halfvec_precache(stats, curr, 0, curr->length-1, curr_R, curr_rd_u, curr_rd_v);
  if(curr_dh_dx == 0.0) return 0.0f;

  const double f_curr = path_measurement_contribution_dx(curr, 0, curr->length-1);
  assert(curr->v[curr->length-1].mode != s_absorb);
  assert(curr->v[curr->length-1].shading.em > 0);

  // rubbish path, sorry dudes.
  if(!(f_curr > 0.0))
  {
    stats->no_throughput++;
    return 0.0f;
  }

  // create a new path by perturbing half vectors:
  if(halfvec_perturb_single(stats, curr, tent, 0, curr->length-1, curr_R, curr_rd_u, curr_rd_v))
  {
    stats->perturb_failed++;
    return 0.0f;
  }

  // TODO: smooth vertices
  double tent_dh_dx = halfvec_precache(stats, tent, 0, tent->length-1, tent_R, tent_rd_u, tent_rd_v);
  if(tent_dh_dx == 0.0) return 0.0f;

  // compute transition probabilities:
  const double p_tent = curr_dh_dx * halfvec_pdf_perturb_single(stats, curr, tent, 0, curr->length-1, 0, curr_R, curr_rd_u, curr_rd_v);
  if(!(p_tent > 0.0f))
  {
    stats->no_throughput++;
    return 0.0f;
  }

  // computing the reverse pdf requires ray differentials to be initialized on tent
  const double p_curr = tent_dh_dx * halfvec_pdf_perturb_single(stats, tent, curr, 0, curr->length-1, 1, tent_R, tent_rd_u, tent_rd_v);
  if(!(p_curr > 0.0f))
  {
    stats->reverse_check_failed++;
    return 0.0f;
  }

  const double f_tent = path_measurement_contribution_dx(tent, 0, tent->length-1);
  if(!(f_tent > 0.0f))
  {
    stats->no_throughput++;
    return 0.0f;
  }
  stats->mutations_proposed++;
  for(int v=1;v<tent->length;v++)
    assert(fabsf(dotproduct(tent->e[v].omega, tent->e[v].omega) - 1.0f) < 1e-4f);

  // compute acceptance as f->/T->  /  f<-/T<-
  return (f_tent/p_tent) / (f_curr/p_curr);
}

void hmc_print_info(FILE *f, void *data)
{
  vmlt_hmc_t *d = (vmlt_hmc_t *)data;
  halfvec_stats_t stats;
  hmc_stats_accum(&stats, d);
  fprintf(f, "         : fake-hmc mutation\n");
  fprintf(f, "           failed project         %.02f%%\n", 100.0f*stats.project_failed/(float)stats.mutations);
  fprintf(f, "           failed perturbation    %.02f%%\n", 100.0f*stats.perturb_failed/(float)stats.mutations);
  fprintf(f, "           no throughput          %.02f%%\n", 100.0f*stats.no_throughput/(float)stats.mutations);
  fprintf(f, "           reverse pdf 0          %.02f%%\n", 100.0f*stats.reverse_check_failed/(float)stats.mutations);
  fprintf(f, "           raydiff compute failed %.02f%%\n", 100.0f*stats.raydiff_compute_failed/(float)stats.mutations);
  fprintf(f, "           successful mutations   %.02f%%\n", 100.0f*stats.mutations_proposed/(float)stats.mutations);
}
