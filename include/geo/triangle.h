#ifndef _CORONA_TRIANGLE_H
#define _CORONA_TRIANGLE_H

#include "corona_common.h"
#include "pathspace.h"

static inline void geo_tri_get_bounds_shutter_open(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  const float *v = geo_get_vertex(p, pi, 0);
  *min = v[dim];
  *max = v[dim];
  for(int k=1;k<vcnt;k++)
  {
    const float *v = geo_get_vertex(p, pi, k);
    *min = fminf(v[dim], *min);
    *max = fmaxf(v[dim], *max);
  }
}

static inline void geo_tri_get_bounds_shutter_close(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  const float *v = geo_get_vertex_shutter_close(p, pi, 0);
  *min = v[dim];
  *max = v[dim];
  for(int k=1;k<vcnt;k++)
  {
    const float *v = geo_get_vertex_shutter_close(p, pi, k);
    *min = fminf(v[dim], *min);
    *max = fmaxf(v[dim], *max);
  }
}

// triangle helper functions for intersection, derivatives, normals, etc.
static inline float geo_tri_get_area(
    const float *v0,
    const float *v1,
    const float *v2)
{
  float e1[3], e2[3], normal[3];
  for(int k=0;k<3;k++)
  {
    e1[k] = v1[k] - v0[k];
    e2[k] = v2[k] - v0[k];
  }
  crossproduct(e1, e2, normal);
  return sqrtf(dotproduct(normal,normal))*.5f;
}

static inline void geo_tri_retime(
    const float *v0,
    const float *v1,
    const float *v2,
    const float u,
    const float v,
    hit_t *hit)
{
  const float w = 1.0f-u-v;
  for(int k=0;k<3;k++) hit->x[k] = w*v0[k] + v*v1[k] + u*v2[k];
}

static inline void geo_tri_get_normal(
    const float *v0,
    const float *v1,
    const float *v2,
    const float *n0,
    const float *n1,
    const float *n2,
    const float u,
    const float v,
    hit_t *hit)
{
  hit->gn[0] = (v1[1] - v0[1])*(v2[2] - v0[2]) - (v1[2] - v0[2])*(v2[1] - v0[1]);
  hit->gn[1] = (v1[2] - v0[2])*(v2[0] - v0[0]) - (v1[0] - v0[0])*(v2[2] - v0[2]);
  hit->gn[2] = (v1[0] - v0[0])*(v2[1] - v0[1]) - (v1[1] - v0[1])*(v2[0] - v0[0]);
  normalise(hit->gn);

  const float w = 1.0f - u - v;
  for(int k=0;k<3;k++) hit->n[k] = u * n2[k] + v * n1[k] + w * n0[k];
  normalise(hit->n);
}

static inline int geo_tri_dpduv(
    const float *v0,
    const float *v1,
    const float *v2,
    const float *uv0,
    const float *uv1,
    const float *uv2,
    const float u,
    const float v,
    float *dpdu,
    float *dpdv)
{
  // compute dpd[uv] wrt the triangle axes:
  for(int k=0;k<3;k++) dpdu[k] = v2[k] - v0[k];
  for(int k=0;k<3;k++) dpdv[k] = v1[k] - v0[k];

  float duv1[2], duv2[2];
  for(int k=0;k<2;k++) duv1[k] = uv2[k] - uv0[k];
  for(int k=0;k<2;k++) duv2[k] = uv1[k] - uv0[k];
  const float det = duv1[0] * duv2[1] - duv2[0] * duv1[1];
  if(det == 0.0f) return 1;
  // transform backwards.
  float dpdu2[3], dpdv2[3];
  for(int k=0;k<3;k++)
    dpdu2[k] = (duv2[1] * dpdu[k] - duv1[1] * dpdv[k])/det;
  for(int k=0;k<3;k++)
    dpdv2[k] = (-duv2[0] * dpdu[k] + duv1[0] * dpdv[k])/det;

  for(int k=0;k<3;k++) dpdu[k] = dpdu2[k];
  for(int k=0;k<3;k++) dpdv[k] = dpdv2[k];
  return 0;
}

static inline void geo_tri_dnduv(
    const float *v0,   // unused
    const float *v1,   // unused
    const float *v2,   // unused
    const float *n0,
    const float *n1,
    const float *n2,
    const float *uv0,
    const float *uv1,
    const float *uv2,
    const float u,
    const float v,
    const float *dpdu,   // unused
    const float *dpdv,   // unused
    float *dndu,
    float *dndv)
{
  const float w = 1.0f - u - v;

  float n[3]; // need length of un-normalised normal, so recompute:
  for(int k=0;k<3;k++) n[k] = u * n2[k] + v * n1[k] + w * n0[k];
  const float invl = 1.f/sqrtf(dotproduct(n, n));
  for(int k=0;k<3;k++) n[k] *= invl;

  // compute dnd[uv] wrt tangent vectors along the triangle axes:
  for(int k=0;k<3;k++) dndu[k] = (n2[k] - n0[k])*invl;
  float dot = dotproduct(n, dndu);
  for(int k=0;k<3;k++) dndu[k] -= n[k] * dot;
  for(int k=0;k<3;k++) dndv[k] = (n1[k] - n0[k])*invl;
  dot = dotproduct(n, dndv);
  for(int k=0;k<3;k++) dndv[k] -= n[k] * dot;

  float duv1[2], duv2[2];
  for(int k=0;k<2;k++) duv1[k] = uv2[k] - uv0[k];
  for(int k=0;k<2;k++) duv2[k] = uv1[k] - uv0[k];
  const float det = duv1[0] * duv2[1] - duv2[0] * duv1[1];
  if(det == 0.0f)
  {
    for(int k=0;k<3;k++) dndu[k] = dndv[k] = 0.0f;
    return;
  }
  // transform backwards.
  float dndu2[3], dndv2[3];
  for(int k=0;k<3;k++)
    dndu2[k] = (duv2[1] * dndu[k] - duv1[1] * dndv[k])/det;
  for(int k=0;k<3;k++)
    dndv2[k] = (-duv2[0] * dndu[k] + duv1[0] * dndv[k])/det;

  for(int k=0;k<3;k++) dndu[k] = dndu2[k];
  for(int k=0;k<3;k++) dndv[k] = dndv2[k];
}



#if 0 // double, conservative chirkov-style
static inline int geo_tri_intersect(
    const float *v0,
    const float *v1,
    const float *v2,
    const primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  long double ac[3], ab[3], bc[3], aO[3], bO[3], cO[3];
  for(int k=0;k<3;k++)
  {
    ab[k] = v1[k] - v0[k];
    ac[k] = v2[k] - v0[k];
    bc[k] = v2[k] - v1[k];
    aO[k] = ray->pos[k] - v0[k];
    bO[k] = ray->pos[k] - v1[k];
    cO[k] = ray->pos[k] - v2[k];
  }
  long double a[3], b[3], c[3];
  // out facing normals for ccw tris:
  crossproduct(ab, aO, a);
  crossproduct(bc, bO, b);
  crossproduct(cO, ac, c);

  // pad 4 double-ulp here, to grow the triangle to avoid holes in between them
  const long double eps = 0.0000000000000004;
  if(((dotproduct(a, ray->dir) >= -eps) && (dotproduct(b, ray->dir) >= -eps) && (dotproduct(c, ray->dir) >= -eps)) ||
     ((dotproduct(a, ray->dir) <= eps) && (dotproduct(b, ray->dir) <= eps) && (dotproduct(c, ray->dir) <= eps)))
  {
    long double n[3];
    crossproduct(ab, ac, n);
    const long double d_inv = 1.0/dotproduct(ray->dir, n);
    const long double dist = -dotproduct(aO, n)*d_inv;
    if(dist > 0 && dist <= hit->dist)
    {
      hit->dist = dist;
      hit->prim = pi;
      hit->u = dotproduct(a, ray->dir)*d_inv;
      hit->v = dotproduct(c, ray->dir)*d_inv;
      return 1;
    }
  }
  return 0;
}
#endif
#if 0 // float, conservative chirkov-style
static inline int geo_tri_intersect(
    const float *v0,
    const float *v1,
    const float *v2,
    const primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  float ac[3], ab[3], bc[3], aO[3], bO[3], cO[3];
  for(int k=0;k<3;k++)
  {
    ab[k] = v1[k] - v0[k];
    ac[k] = v2[k] - v0[k];
    bc[k] = v2[k] - v1[k];
    aO[k] = ray->pos[k] - v0[k];
    bO[k] = ray->pos[k] - v1[k];
    cO[k] = ray->pos[k] - v2[k];
  }
  float a[3], b[3], c[3];
  // out facing normals for ccw tris:
  crossproduct(ab, aO, a);
  crossproduct(bc, bO, b);
  crossproduct(cO, ac, c);

  const float eps = 0;//.00000024f;
  if(((dotproduct(a, ray->dir) >= -eps) && (dotproduct(b, ray->dir) >= -eps) && (dotproduct(c, ray->dir) >= -eps)) ||
     ((dotproduct(a, ray->dir) <= eps) && (dotproduct(b, ray->dir) <= eps) && (dotproduct(c, ray->dir) <= eps)))
  {
    float n[3];
    crossproduct(ab, ac, n);
    const float d_inv = 1.0f/dotproduct(ray->dir, n);
    const float dist = -dotproduct(aO, n)*d_inv;
    if(dist > 0 && dist <= hit->dist)
    {
      hit->dist = dist;
      hit->prim = pi;
      hit->u = dotproduct(a, ray->dir)*d_inv;
      hit->v = dotproduct(c, ray->dir)*d_inv;
      return 1;
    }
  }
  return 0;
}
#endif
#if 1 // float moeller trumbore, epsilon free
static inline int geo_tri_intersect(
    const float *v0,
    const float *v1,
    const float *v2,
    const primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  if(primid_eq(pi, ray->ignore)) return 0;
  float edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  float det,inv_det;

  for(int k=0;k<3;k++)
  {
    edge1[k] = v1[k] - v0[k];
    edge2[k] = v2[k] - v0[k];
  }

  crossproduct(ray->dir, edge2, pvec);
  det = dotproduct(edge1, pvec);
  inv_det = 1.0f / det;

  for(int k=0;k<3;k++) tvec[k] = ray->pos[k] - v0[k];

  float v = dotproduct(tvec, pvec) * inv_det;
  if (v < 0.0f || v > 1.0f) return 0;

  crossproduct(tvec, edge1, qvec);

  float u = dotproduct(ray->dir, qvec) * inv_det;
  if (u < 0.0f || u + v > 1.0f) return 0;

  float dist = dotproduct(edge2, qvec) * inv_det;
  if(dist > ray->min_dist && dist <= hit->dist)
  {
    hit->dist = dist;
    hit->prim = pi;
    hit->u = u;
    hit->v = v;
    return 1;
  }
  return 0;
}
#endif

static inline int geo_tri_intersect_visible(
    const float *v0,
    const float *v1,
    const float *v2,
    const ray_t *ray,
    const float max_dist)
{
  // moller trumbore
  float edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  float det,inv_det;

  for(int k=0;k<3;k++)
  {
    edge1[k] = v1[k] - v0[k];
    edge2[k] = v2[k] - v0[k];
  }

  crossproduct(ray->dir, edge2, pvec);
  det = dotproduct(edge1, pvec);
  inv_det = 1.0 / det;

  for(int k=0;k<3;k++) tvec[k] = ray->pos[k] - v0[k];

  float v = dotproduct(tvec, pvec) * inv_det;
  if (v < 0.0 || v > 1.0) return 0;

  crossproduct(tvec, edge1, qvec);

  float u = dotproduct(ray->dir, qvec) * inv_det;
  if (u < 0.0 || u + v > 1.0) return 0;

  /* calculate t, ray intersects triangle. inv_det calculated and mul'ed for precision reasons :( */
  float dist = dotproduct(edge2, qvec) * inv_det;
  if(dist > 0.0 && dist <= max_dist) return 1;
  return 0;
}



#endif
