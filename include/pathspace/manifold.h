#pragma once

#include "pathspace.h"
#include "matrix2.h"
#include "matrix3.h"
#include "prims.h"
#include "geo.h"
#include "geo/line.h"
#include "pathspace/manifold_vol.h"

static inline float _manifold_inv(const path_t *p, const int v, const float *m, float *inv)
{ // wrapper for fast path for 2x2 blocks
  if(p->v[v].mode & s_volume) return mat3_invert(m, inv);
  else return mat3_invert_sub2(m, inv);
}

// solve for new x positions given h
// returns != 0 on error (singular blocks)
static inline int manifold_map_h_to_x(
    const path_t *curr, // original path to use diffgeo struct and original point position from.
    path_t *tent,       // path under consideration, vertices v[start+1]..v[end-1] will be changed
    const float *dh,    // array of 3*length float representing delta projected half vectors, where length=end-start+1. third component is 0 or distance delta in volumes.
    const float step,   // step size to multiply to vertex offsets
    const int s,        // fixed vertex at start of glossy chain
    const int e)        // fixed vertex at end of glossy chain
{
  // last vertex in an array where [0] means v[s]
  const int n = e-s;

  if(e-s == 2)
  {
    float bi[9], x[3];
    if(_manifold_inv(curr, s+1, curr->v[s+1].diffgeo.b, bi) == 0.0f) return 1;
    mat3_mulv(bi, dh+3, x);
    for(int i=0;i<3;i++) tent->v[s+1].hit.x[i] = curr->v[s+1].hit.x[i] + step * (
      x[0] * curr->v[s+1].diffgeo.dpdu[i] +
      x[1] * curr->v[s+1].diffgeo.dpdv[i]);
    if(tent->v[s+1].mode & s_volume)
      for(int i=0;i<3;i++) tent->v[s+1].hit.x[i] += step * x[2] * curr->v[s+1].hit.gn[i];
    return 0;
  }

  // invert only up to vertex n-1, the last one is fixed and doesn't have a half vector constraint.

  float Li[n+1][9];
  float A[n+1][9];
  float H[n+1][3];
  float tmp[9], tmp2[9];
  if(_manifold_inv(curr, s, curr->v[s].diffgeo.b, Li[0]) == 0.0f) return 1;
  for(int k=1;k<n;k++)
  {
    mat3_mul(curr->v[s+k].diffgeo.a, Li[k-1], A[k]);
    mat3_mul(A[k], curr->v[s+k-1].diffgeo.c, tmp);
    mat3_sub(curr->v[s+k].diffgeo.b, tmp, tmp2);
    // L = b - au
    if(_manifold_inv(curr, s+k, tmp2, Li[k]) == 0.0f) return 1;
  }

  H[0][0] = dh[0]; H[0][1] = dh[1]; H[0][2] = dh[2];
  for(int k=1;k<n;k++)
  {
    mat3_mulv(A[k], H[k-1], tmp);
    for(int i=0;i<3;i++) H[k][i] = dh[3*k+i] - tmp[i];
  }

  float x[n+1][3];
  memset(x, 0, sizeof(float)*(n+1)*3);
  mat3_mulv(Li[n-1], H[n-1], x[n-1]);
  for(int i=0;i<3;i++) tent->v[e-1].hit.x[i] = curr->v[e-1].hit.x[i] + step * (
    x[n-1][0] * curr->v[e-1].diffgeo.dpdu[i] +
    x[n-1][1] * curr->v[e-1].diffgeo.dpdv[i]);
  if(tent->v[e-1].mode & s_volume)
    for(int i=0;i<3;i++) tent->v[e-1].hit.x[i] += step * x[n-1][2] * curr->v[e-1].hit.gn[i];
  for(int k=n-2;k>0;k--)
  {
    mat3_mulv(curr->v[s+k].diffgeo.c, x[k+1], tmp);
    for(int i=0;i<3;i++) tmp[i] = H[k][i] - tmp[i];
    mat3_mulv(Li[k], tmp, x[k]);
    for(int i=0;i<3;i++) tent->v[s+k].hit.x[i] = curr->v[s+k].hit.x[i] + step * (
      x[k][0] * curr->v[s+k].diffgeo.dpdu[i] +
      x[k][1] * curr->v[s+k].diffgeo.dpdv[i]);
    if(tent->v[s+k].mode & s_volume)
      for(int i=0;i<3;i++) tent->v[s+k].hit.x[i] += step * x[k][2] * curr->v[s+k].hit.gn[i];
  }
  return 0;

  // http://faculty.washington.edu/finlayso/ebook/algebraic/advanced/LUtri.htm
  // dk  := h_k (deltas for hu, hv, t)
  // xk  := v[k].x in tangent space
  // b'k := L_k = b_k - a_k*u_{k-1}
  // a'k := a_k * L^{-1}_{k-1}
  //
  // b'1 = b1
  // for k=2..n
  //   a'k = ak/b'{k-1}
  //   b'k = bk - ak/b'{k-1} c{k-1}
  //
  // d'1 = d1
  // for k=2..n
  //   d'k = dk - a'k d'{k-1}
  //
  // xn = d'n / b'n
  // for k=n-1..1
  //    xk = (d'k - ck x{k+1})/b'k
  //
}

// initialize some essential variables in the diffgeo struct on each vertex.
// don't actually compute the tangent space mappings yet.
static inline void manifold_init(path_t *p, int v)
{
  // compute dndu dndv for every vertex, depends on uv coordinates to be useful
  if(!primid_invalid(p->v[v].hit.prim) &&
    (p->v[v].hit.prim.vcnt == 2) &&
     geo_line_is_linestrip(rt.prims, p->v[v].hit.prim))
  { // line strip hair fibers:
    // gn, dpdu, dpdv make up a tangent frame for half vectors, just like for
    // volumes. this means gn = omega[v]
    // n is the hair fiber direction and a and b are not initialised.
    // a,b are unused by the hair scattering bsdf (the frame rotates with
    // incoming direction instead) and in fact should be unused throughout
    float4_t v0, v1;
    geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 0, p->time, &v0);
    geo_get_vertex_time(rt.prims, p->v[v].hit.prim, 1, p->time, &v1);
    for(int k=0;k<3;k++) p->v[v].hit.n[k] = v1.f[k] - v0.f[k];
    normalise(p->v[v].hit.n);
    for(int k=0;k<3;k++) p->v[v].hit.gn[k] = p->e[v].omega[k];
    get_scrambled_onb(p->tangent_frame_scrambling,
        p->v[v].hit.gn, p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv);
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndu[k] = 0.0f;
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] = 0.0f;
    for(int k=0;k<3;k++) p->v[v].hit.a[k] = p->v[v].hit.b[k] = 0.0f/0.0f; // poison
  }
  else if(!primid_invalid(p->v[v].hit.prim))
  { // surface (includes tris, quads, cylinders, spheres)
    // for a significant speedup especially in trivial scenes where
    // such overhead dominates, set this to #if 1:
#ifdef COMPUTE_MANIFOLD_STUFF
    // fill hit.gn and hit.n and hit.{s,t}
    prims_get_normal_time(rt.prims, p->v[v].hit.prim, &p->v[v].hit, p->time);
    // init geo derivatives without flipping the normal (to keep dpdu/dpdv oriented the same way in case
    // the path was bidirectionally constructed the other way around), and do the flipping business
    // manually afterwards (flip normals, dnduv get a sign)
    prims_init_derivatives(rt.prims, p->v[v].hit.prim, &p->v[v].hit, p->time, &p->v[v].diffgeo);

    // flip shading normals towards ray and set inside bit
    if(dotproduct(p->e[v].omega, p->v[v].hit.gn) > 0.0f)
    {
      for(int i=0;i<3;i++) p->v[v].hit.n[i] = - p->v[v].hit.n[i];
      p->v[v].flags |= s_inside;
      for(int k=0;k<3;k++) p->v[v].diffgeo.dndu[k] *= -1.0f;
      for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] *= -1.0f;
    }
    else p->v[v].flags &= ~s_inside;

    // ortho normalise, stick to dpdu as direction:
    const float ilu = 1.0f/sqrtf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdu));
    for(int k=0;k<3;k++) p->v[v].diffgeo.dpdu[k] *= ilu;
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndu[k] *= ilu;
    const float dotu = dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv);
    for(int k=0;k<3;k++) p->v[v].diffgeo.dpdv[k] -= dotu*p->v[v].diffgeo.dpdu[k];
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] -= dotu*p->v[v].diffgeo.dndu[k];
    const float dotn = dotproduct(p->v[v].hit.gn, p->v[v].diffgeo.dpdv);
    for(int k=0;k<3;k++) p->v[v].diffgeo.dpdv[k] -= dotn*p->v[v].hit.gn[k];
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] -= dotn*p->v[v].hit.gn[k];
    const float ilv = 1.0f/sqrtf(dotproduct(p->v[v].diffgeo.dpdv, p->v[v].diffgeo.dpdv));
    for(int k=0;k<3;k++) p->v[v].diffgeo.dpdv[k] *= ilv;
    for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] *= ilv;

    // damnit. this happens for super thin anisotropic quads:
    if(fabsf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv)) >= 1e-4f)
      get_onb(p->v[v].hit.gn, p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv);
    assert(fabsf(dotproduct(p->v[v].hit.gn,       p->v[v].hit.gn) - 1.0f) < 1e-4f);
    if(fabsf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdu) - 1.0f) >= 1e-4f)
    {
      fprintf(stderr, "normal %g %g %g du %g %g %g dv %g %g %g\n",
          p->v[v].hit.gn[0],
          p->v[v].hit.gn[1],
          p->v[v].hit.gn[2],
          p->v[v].diffgeo.dpdu[0],
          p->v[v].diffgeo.dpdu[1],
          p->v[v].diffgeo.dpdu[2],
          p->v[v].diffgeo.dpdv[0],
          p->v[v].diffgeo.dpdv[1],
          p->v[v].diffgeo.dpdv[2]);
    }
    assert(fabsf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdu) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].diffgeo.dpdv, p->v[v].diffgeo.dpdv) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].diffgeo.dpdv)) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].diffgeo.dpdu, p->v[v].hit.gn)) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].diffgeo.dpdv, p->v[v].hit.gn)) < 1e-4f);

    // shading frame is also rotated to match dpdu:
    float dot_dpdu_n = dotproduct(p->v[v].diffgeo.dpdu, p->v[v].hit.n);
    for(int k=0;k<3;k++) p->v[v].hit.a[k] = p->v[v].diffgeo.dpdu[k] - dot_dpdu_n * p->v[v].hit.n[k];
    normalise(p->v[v].hit.a);

    float dot_dpdv_n = dotproduct(p->v[v].diffgeo.dpdv, p->v[v].hit.n);
    for(int k=0;k<3;k++) p->v[v].hit.b[k] = p->v[v].diffgeo.dpdv[k] - dot_dpdv_n * p->v[v].hit.n[k];
    float dot_a_b = dotproduct(p->v[v].hit.a, p->v[v].hit.b);
    for(int k=0;k<3;k++) p->v[v].hit.b[k] -= dot_a_b * p->v[v].hit.a[k];
    normalise(p->v[v].hit.b);
    // unfortunately sometimes dot_dpdv_n ~= 1.0 because <n,gn>~=0 and then this happens:
    if(fabsf(dotproduct(p->v[v].hit.a, p->v[v].hit.b)) >= 1e-4f || 
       fabsf(dotproduct(p->v[v].hit.b, p->v[v].hit.n)) >= 1e-4f)
      get_onb(p->v[v].hit.n, p->v[v].hit.a, p->v[v].hit.b);

    assert(fabsf(dotproduct(p->v[v].hit.a, p->v[v].hit.a) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].hit.b, p->v[v].hit.b) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].hit.n, p->v[v].hit.n) - 1.0f) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].hit.a, p->v[v].hit.b)) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].hit.a, p->v[v].hit.n)) < 1e-4f);
    assert(fabsf(dotproduct(p->v[v].hit.b, p->v[v].hit.n)) < 1e-4f);
#else
    prims_get_normal_time(rt.prims, p->v[v].hit.prim, &p->v[v].hit, p->time);
    // flip shading normals towards ray and set inside bit
    if(dotproduct(p->e[v].omega, p->v[v].hit.gn) > 0.0f)
    {
      for(int i=0;i<3;i++) p->v[v].hit.n[i] = - p->v[v].hit.n[i];
      p->v[v].flags |= s_inside;
      for(int k=0;k<3;k++) p->v[v].diffgeo.dndu[k] *= -1.0f;
      for(int k=0;k<3;k++) p->v[v].diffgeo.dndv[k] *= -1.0f;
    }
    else p->v[v].flags &= ~s_inside;
    get_scrambled_onb(p->tangent_frame_scrambling, p->v[v].hit.n, p->v[v].hit.a, p->v[v].hit.b);
    for(int k=0;k<3;k++)
    {
      p->v[v].diffgeo.dpdu[k] = p->v[v].hit.a[k];
      p->v[v].diffgeo.dpdv[k] = p->v[v].hit.b[k];
      p->v[v].diffgeo.dndu[k] = 0.0f;
      p->v[v].diffgeo.dndv[k] = 0.0f;
    }
#endif
  }
  else
  { // volume scattering
    for(int k=0;k<3;k++) p->v[v].hit.n[k] = p->v[v].hit.gn[k] = p->e[v].omega[k];
    get_scrambled_onb(p->tangent_frame_scrambling, p->v[v].hit.n, p->v[v].hit.a, p->v[v].hit.b);
    for(int k=0;k<3;k++)
    {
      p->v[v].diffgeo.dpdu[k] = p->v[v].hit.a[k];
      p->v[v].diffgeo.dpdv[k] = p->v[v].hit.b[k];
      p->v[v].diffgeo.dndu[k] = 0.0f;
      p->v[v].diffgeo.dndv[k] = 0.0f;
    }
  }
}

// compute the half vector derivative matrices a, b, c for the given vertex p->v[i].
static inline int manifold_compute_derivatives(path_t *p, int i)
{
  vertex_manifold_t *mf = &p->v[i].diffgeo;
  float wo[3]; // pointing from vert i to i+1
  for(int k=0;k<3;k++) wo[k] = p->v[i+1].hit.x[k] - p->v[i].hit.x[k];
  float ilo = sqrtf(dotproduct(wo, wo));
  if(ilo == 0.0f) return 1;
  ilo = 1.0/ilo;
  for(int k=0;k<3;k++) wo[k] *= ilo;

  // init blocks outside 2x2 sub block to 0
  mat3_set_zero(mf->a);
  mat3_set_identity(mf->b);
  mat3_set_zero(mf->c);

  if(mf->type == s_pinned_position)
  {
    return 0;
  }
  else if(mf->type == s_pinned_direction)
  {
    float dC_dnext_u[3], dC_dnext_v[3], dC_dcurr_u[3], dC_dcurr_v[3];
    const float dot_wo_v0dpdu = dotproduct(wo, p->v[i].diffgeo.dpdu);
    const float dot_wo_v0dpdv = dotproduct(wo, p->v[i].diffgeo.dpdv);
    const float dot_wo_v1dpdu = dotproduct(wo, p->v[i+1].diffgeo.dpdu);
    const float dot_wo_v1dpdv = dotproduct(wo, p->v[i+1].diffgeo.dpdv);
    for(int k=0;k<3;k++)
    {
      dC_dnext_u[k] = (p->v[i+1].diffgeo.dpdu[k] - wo[k] * dot_wo_v1dpdu) * ilo;
      dC_dnext_v[k] = (p->v[i+1].diffgeo.dpdv[k] - wo[k] * dot_wo_v1dpdv) * ilo;
      dC_dcurr_u[k] = (wo[k] * dot_wo_v0dpdu - p->v[i].diffgeo.dpdu[k]) * ilo;
      dC_dcurr_v[k] = (wo[k] * dot_wo_v0dpdv - p->v[i].diffgeo.dpdv[k]) * ilo;
    }

    mf->b[0] = dotproduct(dC_dcurr_u, p->v[0].diffgeo.dpdu);
    mf->b[1] = dotproduct(dC_dcurr_v, p->v[0].diffgeo.dpdu);
    mf->b[3] = dotproduct(dC_dcurr_u, p->v[0].diffgeo.dpdv);
    mf->b[4] = dotproduct(dC_dcurr_v, p->v[0].diffgeo.dpdv);
    mf->b[8] = 1.0f;

    mf->c[0] = dotproduct(dC_dnext_u, p->v[0].diffgeo.dpdu);
    mf->c[1] = dotproduct(dC_dnext_v, p->v[0].diffgeo.dpdu);
    mf->c[3] = dotproduct(dC_dnext_u, p->v[0].diffgeo.dpdv);
    mf->c[4] = dotproduct(dC_dnext_v, p->v[0].diffgeo.dpdv);
    return 0;
  }

  float wi[3]; // pointing from vert i to i-1
  for(int k=0;k<3;k++) wi[k] = p->v[i-1].hit.x[k] - p->v[i].hit.x[k];
  float ili = sqrtf(dotproduct(wi, wi));
  if(ili == 0.0f) return 1;
  ili = 1.0/ili;
  for(int k=0;k<3;k++) wi[k] *= ili;

  if(p->v[i].mode & (s_reflect | s_transmit)) /* reflection or refraction */
  {
    // ratio n1/n2 where n2 is inside the shape
    const float eta = (p->v[i].mode & s_reflect) ? 1.0f : (1.0f/mf(path_eta_ratio(p, i), 0));
    if(!(eta > 0.0)) return 1; // broken volume nesting
    // reflected path also goes through normalise_H, but happens to have eta == 1.
    const int normalise_H = !((p->v[i].mode & s_transmit) && (eta == 1.0f));

    float ilh, H[3];
    if(normalise_H)
    {
      for(int k=0;k<3;k++) H[k] = wi[k] + eta * wo[k];
      ilh = 1.0f/dotproduct(H, p->v[i].hit.n);
      for(int k=0;k<3;k++) H[k] *= ilh;
    }
    else
    {
      // don't normalise H for index-matched transmission
      for(int k=0;k<3;k++) H[k] = wi[k] + wo[k];
      ilh = 1.0f;
    }

    float dot_H_n    = dotproduct(p->v[i].hit.n, H);
    float dot_H_dndu = dotproduct(mf->dndu, H);
    float dot_H_dndv = dotproduct(mf->dndv, H);
    float dot_dpdu_n = dotproduct(mf->dpdu, p->v[i].hit.n);
    float dot_dpdv_n = dotproduct(mf->dpdv, p->v[i].hit.n);

    // tangent frame around shading normal:
    float *s = p->v[i].hit.a, *t = p->v[i].hit.b;

    ilo *= eta * ilh; ili *= ilh;

    // derivatives of C wrt previous vertex v[i-1]
    float dH_du[3], dH_dv[3], dH_dn[3];
    const float dot_wi_dpdu = dotproduct(wi, p->v[i-1].diffgeo.dpdu);
    const float dot_wi_dpdv = dotproduct(wi, p->v[i-1].diffgeo.dpdv);
    for(int k=0;k<3;k++) dH_du[k] = (p->v[i-1].diffgeo.dpdu[k] - wi[k] * dot_wi_dpdu) * ili;
    for(int k=0;k<3;k++) dH_dv[k] = (p->v[i-1].diffgeo.dpdv[k] - wi[k] * dot_wi_dpdv) * ili;

    if(normalise_H)
    {
      const float dotu = dotproduct(dH_du, p->v[i].hit.n);
      const float dotv = dotproduct(dH_dv, p->v[i].hit.n);
      for(int k=0;k<3;k++)
      {
        dH_du[k] -= H[k] * dotu;
        dH_dv[k] -= H[k] * dotv;
      }
    }

    mf->a[0] = dotproduct(dH_du, s);
    mf->a[1] = dotproduct(dH_dv, s);
    mf->a[3] = dotproduct(dH_du, t);
    mf->a[4] = dotproduct(dH_dv, t);

    if(p->v[i-1].mode & s_volume)
    { // 3x2 case
      const float dot_wi_dpdn = dotproduct(wi, p->v[i-1].hit.gn);
      for(int k=0;k<3;k++) dH_dn[k] = (p->v[i-1].hit.gn[k] - wi[k] * dot_wi_dpdn) * ili;
      if(normalise_H)
      {
        const float dotn = dotproduct(dH_dn, p->v[i].hit.n);
        for(int k=0;k<3;k++) dH_dn[k] -= H[k] * dotn;
      }
      mf->a[2] = dotproduct(dH_dn, s);
      mf->a[5] = dotproduct(dH_dn, t);
    }

    // derivatives of C wrt current vertex v[i]
    for(int k=0;k<3;k++)
    {
      dH_du[k] = -p->v[i].diffgeo.dpdu[k] * (ili + ilo) + wi[k] * dotproduct(wi, p->v[i].diffgeo.dpdu) * ili
        + wo[k] * dotproduct(wo, p->v[i].diffgeo.dpdu) * ilo;
      dH_dv[k] = -p->v[i].diffgeo.dpdv[k] * (ili + ilo) + wi[k] * dotproduct(wi, p->v[i].diffgeo.dpdv) * ili
        + wo[k] * dotproduct(wo, p->v[i].diffgeo.dpdv) * ilo;
    }
    if(normalise_H)
    {
      const float dotu = dot_H_dndu + dotproduct(dH_du, p->v[i].hit.n);
      const float dotv = dot_H_dndv + dotproduct(dH_dv, p->v[i].hit.n);
      for(int k=0;k<3;k++)
      {
        dH_du[k] -= H[k] * dotu;
        dH_dv[k] -= H[k] * dotv;
      }
    }

    mf->b[0] = dotproduct(dH_du, s) - dotproduct(mf->dpdu, mf->dndu) * dot_H_n - dot_dpdu_n * dot_H_dndu;
    mf->b[1] = dotproduct(dH_dv, s) - dotproduct(mf->dpdu, mf->dndv) * dot_H_n - dot_dpdu_n * dot_H_dndv;
    mf->b[3] = dotproduct(dH_du, t) - dotproduct(mf->dpdv, mf->dndu) * dot_H_n - dot_dpdv_n * dot_H_dndu;
    mf->b[4] = dotproduct(dH_dv, t) - dotproduct(mf->dpdv, mf->dndv) * dot_H_n - dot_dpdv_n * dot_H_dndv;
    mf->b[8] = 1.0f;

    // derivatives of C wrt v[i+1]
    for(int k=0;k<3;k++)
    {
      dH_du[k] = (p->v[i+1].diffgeo.dpdu[k] - wo[k] * dotproduct(wo, p->v[i+1].diffgeo.dpdu)) * ilo;
      dH_dv[k] = (p->v[i+1].diffgeo.dpdv[k] - wo[k] * dotproduct(wo, p->v[i+1].diffgeo.dpdv)) * ilo;
    }

    if(normalise_H)
    {
      const float dotu = dotproduct(dH_du, p->v[i].hit.n);
      const float dotv = dotproduct(dH_dv, p->v[i].hit.n);
      for(int k=0;k<3;k++)
      {
        dH_du[k] -= H[k] * dotu;
        dH_dv[k] -= H[k] * dotv;
      }
    }

    // c[2] and c[5] are always zero in this tangent frame at [i+1],
    // even if the next vertex is a 3d volume vertex.
    mf->c[0] = dotproduct(dH_du, s);
    mf->c[1] = dotproduct(dH_dv, s);
    mf->c[3] = dotproduct(dH_du, t);
    mf->c[4] = dotproduct(dH_dv, t);

    // store microfacet normal in plane/plane
    mf->h[0] = dotproduct(H, p->v[i].hit.a);
    mf->h[1] = dotproduct(H, p->v[i].hit.b);
    mf->h[2] = 0.0f; // dummy constraint
  }
  else if(p->v[i].mode & (s_volume | s_fiber))
  { // medium scattering. fibers will only get a similar half vector + tangent frame.
    // 3x3 matrices for these vertices, as they come with an additional distance constraint
    //
    //   xu xv xn
    //  | .  .   | hu
    //  | .  .   | hv
    //  |        | t
    //
    // t is sum of distance k,k-1 and k,k+1
    for(int k=0;k<3;k++) wi[k] = p->v[i].hit.x[k] - p->v[i-1].hit.x[k];
    float lwi = sqrtf(dotproduct(wi, wi));
    for(int k=0;k<3;k++) wo[k] = p->v[i+1].hit.x[k] - p->v[i].hit.x[k];
    float lwo = sqrtf(dotproduct(wo, wo));

    mf->a[0] = manifold_vol_dhu_du('a', wi, wo, lwi, lwo, p->v[i-1].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->a[1] = manifold_vol_dhu_du('a', wi, wo, lwi, lwo, p->v[i-1].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->a[3] = manifold_vol_dhv_du('a', wi, wo, lwi, lwo, p->v[i-1].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->a[4] = manifold_vol_dhv_du('a', wi, wo, lwi, lwo, p->v[i-1].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->a[6] = - dotproduct(p->v[i-1].diffgeo.dpdu, wi)/lwi;
    mf->a[7] = - dotproduct(p->v[i-1].diffgeo.dpdv, wi)/lwi;
    if(p->v[i-1].mode & s_volume)
    { // else a 2x3 block
      mf->a[2] = manifold_vol_dhu_du('a', wi, wo, lwi, lwo, p->v[i-1].hit.gn, p->tangent_frame_scrambling);
      mf->a[5] = manifold_vol_dhv_du('a', wi, wo, lwi, lwo, p->v[i-1].hit.gn, p->tangent_frame_scrambling);
      mf->a[8] = - dotproduct(p->v[i-1].hit.gn, wi)/lwi;
    }

    mf->b[0] = manifold_vol_dhu_du('b', wi, wo, lwi, lwo, p->v[i].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->b[1] = manifold_vol_dhu_du('b', wi, wo, lwi, lwo, p->v[i].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->b[2] = manifold_vol_dhu_du('b', wi, wo, lwi, lwo, p->v[i].hit.n, p->tangent_frame_scrambling);
    mf->b[3] = manifold_vol_dhv_du('b', wi, wo, lwi, lwo, p->v[i].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->b[4] = manifold_vol_dhv_du('b', wi, wo, lwi, lwo, p->v[i].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->b[5] = manifold_vol_dhv_du('b', wi, wo, lwi, lwo, p->v[i].hit.n, p->tangent_frame_scrambling);
    mf->b[6] = 0;// - dotproduct(p->v[i].diffgeo.dpdu, wo)/lwo;
    mf->b[7] = 0;// - dotproduct(p->v[i].diffgeo.dpdv, wo)/lwo;
    mf->b[8] = 1;// - dotproduct(p->v[i].hit.gn, wo)/lwo;

    mf->c[0] = manifold_vol_dhu_du('c', wi, wo, lwi, lwo, p->v[i+1].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->c[1] = manifold_vol_dhu_du('c', wi, wo, lwi, lwo, p->v[i+1].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->c[3] = manifold_vol_dhv_du('c', wi, wo, lwi, lwo, p->v[i+1].diffgeo.dpdu, p->tangent_frame_scrambling);
    mf->c[4] = manifold_vol_dhv_du('c', wi, wo, lwi, lwo, p->v[i+1].diffgeo.dpdv, p->tangent_frame_scrambling);
    mf->c[6] = 0;
    mf->c[7] = 0;
    if(p->v[i+1].mode & s_volume)
    { // else a 2x3 block
      mf->c[2] = 0.0f;//manifold_vol_dhu_du('c', wi, wo, lwi, lwo, p->v[i+1].hit.gn, p->tangent_frame_scrambling);
      mf->c[5] = 0.0f;//manifold_vol_dhv_du('c', wi, wo, lwi, lwo, p->v[i+1].hit.gn, p->tangent_frame_scrambling);
      mf->c[8] = 0.0f;
    }

    // hit->gn == hit->n == wi/|wi|
    mf->h[0] = dotproduct(p->v[i].diffgeo.dpdu, wo)/lwo;
    mf->h[1] = dotproduct(p->v[i].diffgeo.dpdv, wo)/lwo;
    mf->h[2] = dotproduct(p->v[i].hit.gn, wo)/lwo;
    // riemann sphere:
    if(mf->h[2] < -0.99999f) mf->h[0] = mf->h[1] = 0.0f;
    else
    {
      mf->h[0] /= (mf->h[2] + 1.0f);
      mf->h[1] /= (mf->h[2] + 1.0f);
    }
    mf->h[2] = lwi;// + lwo; // don't use edge.dist to only depend on vertex locations
  }
  return 0;
}

// similar to wenzel's code in mitsuba.
// initialize tangent space derivatives dh/dx
// start: diffuse vertex before the spec chain
// end: diffuse vertex at end of spec chain
// depends on valid vertex mode, diffgeo.d[pn]d[uv], hit.x, hit.n (and =hit.gn for volumes) and consistent volume nesting (path_eta_ratio).
// writes diffgeo.{a,b,c,h}
static inline double manifold_compute_tangents(path_t *p, int start, int end)
{
  assert(start >= 0);
  assert(end <= p->length-1);

  mat2_set_zero(p->v[start].diffgeo.Tp);   // does not move at all
  mat2_set_identity(p->v[end].diffgeo.Tp); // moves freely along its own tangent space

  if(end-start == 1) return 1.0f; // all good.

  // assemble matrices a, b, c (how does vertex i depend on i-1, i, and i+1)
  for(int i=start;i<end;i++)
  {
    if(manifold_compute_derivatives(p, i))
      return 0.0;
  }
  return 1.0;
}

