#ifndef CORONA_GEO_LINE_H
#define CORONA_GEO_LINE_H

#include "corona_common.h"
#include "pathspace.h"

// line segments are really truncated cones with two radii r0 and r1 at vertices v0 and v1.

static inline void geo_line_get_radii(const prims_t *p, primid_t pi, float *r0, float *r1)
{
  // could motion blur this, currently we're just using the shutter open radius:
  prims_shape_t *s = p->shape + pi.shapeid;
  *r0 = tofloat(s->vtx[(pi.mb+1)*s->vtxidx[pi.vi + 0].v].n);
  *r1 = tofloat(s->vtx[(pi.mb+1)*s->vtxidx[pi.vi + 1].v].n);
}

static inline int geo_line_is_linestrip(const prims_t *p, primid_t pi)
{
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  return MAX(r0,r1) <= 1e-2f;
}

static inline void geo_line_get_bounds_shutter_open(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const float *v0 = geo_get_vertex(p, pi, 0);
  const float *v1 = geo_get_vertex(p, pi, 1);
  float d[3], a[3], b[3];
  for(int k=0;k<3;k++) d[k] = v1[k] - v0[k];
  normalise(d);
  get_onb(d, a, b);
  const float theta = atan2f(a[dim], b[dim]);
  const float m = fabsf(sinf(theta)*a[dim]) + fabsf(cosf(theta)*b[dim]);
  *min = fminf(v0[dim] - r0*m, v1[dim] - r1*m);
  *max = fmaxf(v0[dim] + r0*m, v1[dim] + r1*m);
}

static inline void geo_line_get_bounds_shutter_close(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const float *v0 = geo_get_vertex_shutter_close(p, pi, 0);
  const float *v1 = geo_get_vertex_shutter_close(p, pi, 1);
  float d[3], a[3], b[3];
  for(int k=0;k<3;k++) d[k] = v1[k] - v0[k];
  normalise(d);
  get_onb(d, a, b);
  const float theta = atan2f(a[dim], b[dim]);
  const float m = fabsf(sinf(theta)*a[dim]) + fabsf(cosf(theta)*b[dim]);
  *min = fminf(v0[dim] - r0*m, v1[dim] - r1*m);
  *max = fmaxf(v0[dim] + r0*m, v1[dim] + r1*m);
}

static inline float geo_line_get_area(const prims_t *p, primid_t pi)
{
  // area is at shutter open
  const float *v0 = geo_get_vertex(p, pi, 0);
  const float *v1 = geo_get_vertex(p, pi, 1);
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const float d[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
  const float h = sqrtf(dotproduct(d, d));
  const float l = sqrtf(r0*r0 + h*h);
  const float a0 = M_PI * r0 * l;
  const float a1 = M_PI * r1 * l;
  return a1 - a0;
}

static inline float geo_line_get_area_time(const prims_t *p, primid_t pi, const float time)
{
  // area is at shutter open
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, time, &v0);
  geo_get_vertex_time(p, pi, 1, time, &v1);
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const float d[3] = {v1.f[0]-v0.f[0], v1.f[1]-v0.f[1], v1.f[2]-v0.f[2]};
  const float h = sqrtf(dotproduct(d, d));
  const float l = sqrtf(r0*r0 + h*h);
  const float a0 = M_PI * r0 * l;
  const float a1 = M_PI * r1 * l;
  return a1 - a0;
}

static inline void geo_line_retime(
    const prims_t *p,
    primid_t pi,
    hit_t *hit,
    const float time)
{
  // sample a point on the cone. u is along the hair, v is around it.
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, time, &v0);
  geo_get_vertex_time(p, pi, 1, time, &v1);
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);

  // we isotropically sample circles at height y, where
  // p(y) = (r0 + (r1-r0)y) * 2/(r1+r0)
  float y;
  if(fabsf(r1-r0) < 1e-3f)
    y = hit->u;
  else
    y = (sqrtf((r1*r1 - r0*r0)*hit->u + r0*r0) - r0)/(r1-r0);
  // sample angle phi around the trunc cone
  const float phi = 2.0*M_PI*hit->v;
  float sinphi, cosphi;
  common_sincosf(phi, &sinphi, &cosphi);

  float d[3], a[3], b[3];
  for(int k=0;k<3;k++) d[k] = v1.f[k] - v0.f[k];
  const float ilen_d = 1.0f/sqrtf(dotproduct(d, d));
  for(int k=0;k<3;k++) d[k] *= ilen_d;
  get_onb(d, a, b);

  // compute vertex:
  for(int k=0;k<3;k++)
    hit->x[k] = v0.f[k] + (v1.f[k] - v0.f[k])*y + a[k] * sinphi + b[k] * cosphi;
}

static inline void geo_line_get_normal_time(const prims_t *p, primid_t pi, hit_t *hit, float time)
{
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, time, &v0);
  geo_get_vertex_time(p, pi, 1, time, &v1);
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  if(fabsf(r0-r1) < 1e-3f && r0 < 0.01f)
  { // hair fibers can't resolve normals at floating point precision, set to 0
    for(int k=0;k<3;k++) hit->n[k] = hit->gn[k] = 0.0f;
    return;
  }
  float d[3], a[3], b[3];
  for(int k=0;k<3;k++) d[k] = v1.f[k] - v0.f[k];
  const float ilen_d = 1.0f/sqrtf(dotproduct(d, d));
  for(int k=0;k<3;k++) d[k] *= ilen_d;
  get_onb(d, a, b);
  const float phi = 2.0*M_PI*hit->v;
  float sinphi, cosphi;
  common_sincosf(phi, &sinphi, &cosphi);
  float n[3];
  for(int k=0;k<3;k++)
    n[k] = a[k] * sinphi + b[k] * cosphi;

  // compute hit normal from u/v (actually only depends on v, u will only move the point on the tangent plane then).
  const float rr = r1-r0;
  if(fabsf(rr) < 1e-3)
  {
    for(int k=0;k<3;k++) hit->n[k] = n[k];
    memcpy(hit->gn, hit->n, sizeof(float)*3);
    return;
  }

  // tilt to match cylinder
  for(int k=0;k<3;k++)
    hit->n[k] = n[k] - d[k] * (r1-r0)*ilen_d;
  normalise(hit->n);
  memcpy(hit->gn, hit->n, sizeof(float)*3);
}

static inline void geo_line_init_derivatives(
    const prims_t *p,
    const primid_t pi,
    const hit_t *hit,
    const float time,
    vertex_manifold_t *mf)
{
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, time, &v0);
  geo_get_vertex_time(p, pi, 1, time, &v1);
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  float d[3], a[3], b[3];
  for(int k=0;k<3;k++) d[k] = v1.f[k] - v0.f[k];
  const float len_d = sqrtf(dotproduct(d, d));
  for(int k=0;k<3;k++) d[k] /= len_d;
  get_onb(d, a, b);

  // u is along the hair fiber, v is around it (0..1 remapped from 0..2pi)
  // u: v0..v1, v: b->a
  const float r = r0 + hit->u*(r1-r0);

  float dpdu[3], dpdv[3];
  // could use hit->x projected to a,b,d instead of sincosf
  float cos_phi, sin_phi;
  common_sincosf(2.0f*M_PI*hit->v, &sin_phi, &cos_phi);
  // clamp so that tip of cone still has some non-zero dpdv
  dpdv[0] = -sin_phi * fmaxf(1e-8f, r) * 2.0f*M_PI;
  dpdv[1] =  cos_phi * fmaxf(1e-8f, r) * 2.0f*M_PI;
  dpdv[2] =  0.0f;
  dpdu[0] = (r1-r0) * cos_phi;
  dpdu[1] = (r1-r0) * sin_phi;
  dpdu[2] = len_d;

  // transform back into world space:
  for(int k=0;k<3;k++)
  {
    mf->dpdu[k] = dpdu[0] * b[k] + dpdu[1] * a[k] + dpdu[2] * d[k];
    mf->dpdv[k] = dpdv[0] * b[k] + dpdv[1] * a[k] + dpdv[2] * d[k];
  }

  for(int k=0;k<3;k++)
  {
    mf->dndv[k] = mf->dpdv[k] / fmaxf(1e-8f, r);
    mf->dndu[k] = 0.0f;
  }
}

#if 0
static inline float _geo_line_intersect_cylinder(
    const float *v0,
    const float *v1,
    const float r,
    const ray_t *ray,
    float *out,
    float *len)
{
  double d[3], a[3], b[3], o[3] = {0.0};
  double dir[3] = {ray->dir[0], ray->dir[1], ray->dir[2]};
  for(int k=0;k<3;k++) d[k] = (double)v1[k] - v0[k];

  if(r < 0.01)
  { // hair intersection code is just line strip (cylinder not precise enough)
    crossproduct(d, dir, a);
    for(int k=0;k<3;k++) o[k] = v0[k] - (double)ray->pos[k];
    const double dotp = dotproduct(o, a);
    const double ilen = 1.0/sqrt(dotproduct(a, a));
    const double dist = fabs(dotp*ilen);
    if(dist > (double)r) return -1.0f;

    // also do line stripe intersection:
    const double dlen = sqrt(dotproduct(d,d));
    if(len) *len = dlen;
    crossproduct(d, a, b);

    const double t = dotproduct(o, b)/dotproduct(b, dir);
    for(int k=0;k<3;k++) out[k] = t * dir[k] - o[k];
    out[0] = dotproduct(out, d)/dlen;
    out[1] = 0.0f;
    out[2] = 1.0f;
    if(out[0] >= 0.0 && out[0] <= dlen)
      return t;
    return -1.0f;
  }

  const double dlen = sqrt(dotproduct(d,d));
  if(len) *len = dlen;
  for(int k=0;k<3;k++) d[k] *= 1.0/dlen;
  if(fabs(d[1]) < 0.5)
  {
    double up[3] = {0, 1, 0};
    crossproduct(d, up, a);
  }
  else
  {
    double rg[3] = {1, 0, 0};
    crossproduct(d, rg, a);
  }
  crossproduct(d, a, b);
  double ilena = 1./sqrt(dotproduct(a, a)), ilenb = 1./sqrt(dotproduct(b,b));
  for(int k=0;k<3;k++) a[k] *= ilena;
  for(int k=0;k<3;k++) b[k] *= ilenb;
  // convert ray to cylinder space, where x-axis is aligned with d:
  double w[3] = {0.0};
  for(int k=0;k<3;k++)
  {
    o[0] += ((double)ray->pos[k] - v0[k])*d[k];
    o[1] += ((double)ray->pos[k] - v0[k])*a[k];
    o[2] += ((double)ray->pos[k] - v0[k])*b[k];
    w[0] += dir[k]*d[k];
    w[1] += dir[k]*a[k];
    w[2] += dir[k]*b[k];
  }
  // intersect ray with circle around yz:
  const double A = w[1]*w[1] + w[2]*w[2];
  const double B = 2.0f*(o[1]*w[1]+o[2]*w[2]);
  const double C = o[1]*o[1]+o[2]*o[2] - r*r;
  const double discr = B*B - 4.0*A*C;
  if(discr < 0.0) return -1.0f;
  double temp, sqrt_discrim = sqrt(discr);
  if (B < 0)
    temp = -0.5 * (B - sqrt_discrim);
  else
    temp = -0.5 * (B + sqrt_discrim);
  const double t0 = temp / A;
  const double t1 = C / temp;

  double t;
  if(t0 <= 0.0) t = t1;
  else if(t1 <= 0.0) t = t0;
  else
  {
    t = fmin(t0, t1);
    for(int i=0;i<2;i++)
    {
      for(int k=0;k<3;k++) out[k] = o[k] + t*w[k];
      if(out[0] >= 0.0 && out[0] <= dlen)
        return t;
      t = fmax(t0, t1);
    }
    return -1.0f;
  }

  for(int k=0;k<3;k++)
    out[k] = o[k] + t*w[k];
  if(out[0] >= 0.0 && out[0] <= dlen)
    return t;
  else return -1.0f;
}
#else
static inline float _geo_line_intersect_cylinder(
    const float *v0,
    const float *v1,
    const float r,
    const ray_t *ray,
    float *out,
    float *len)
{
  float d[3], a[3], b[3], o[3] = {0.0};
  for(int k=0;k<3;k++) d[k] = v1[k] - v0[k];

  if(r < 0.01)
  { // hair intersection code is just line strip (cylinder not precise enough)
    crossproduct(d, ray->dir, a);
    for(int k=0;k<3;k++) o[k] = v0[k] - ray->pos[k];
    const float dotp = dotproduct(o, a);
    const float ilen = 1.0/sqrtf(dotproduct(a, a));
    const float dist = fabsf(dotp*ilen);
    if(dist > r) return -1.0f;

    // also do line stripe intersection:
    const float dlen = sqrtf(dotproduct(d,d));
    if(len) *len = dlen;
    crossproduct(d, a, b);

    const float t = dotproduct(o, b)/dotproduct(b, ray->dir);
    for(int k=0;k<3;k++) out[k] = t * ray->dir[k] - o[k];
    out[0] = dotproduct(out, d)/dlen;
    out[1] = 0.0f;
    out[2] = 1.0f;
    if(out[0] >= 0.0 && out[0] <= dlen)
      return t;
    return -1.0f;
  }

  const float dlen = sqrtf(dotproduct(d,d));
  if(len) *len = dlen;
  for(int k=0;k<3;k++) d[k] *= 1.0f/dlen;
  get_onb(d, a, b);
  // convert ray to cylinder space, where x-axis is aligned with d:
  float w[3] = {0.0};
  for(int k=0;k<3;k++)
  {
    o[0] += (ray->pos[k] - v0[k])*d[k];
    o[1] += (ray->pos[k] - v0[k])*a[k];
    o[2] += (ray->pos[k] - v0[k])*b[k];
    w[0] += ray->dir[k]*d[k];
    w[1] += ray->dir[k]*a[k];
    w[2] += ray->dir[k]*b[k];
  }
  // intersect ray with circle around yz:
  const float A = w[1]*w[1] + w[2]*w[2];
  const float B = 2.0f*(o[1]*w[1]+o[2]*w[2]);
  const float C = o[1]*o[1]+o[2]*o[2] - r*r;
  const float discr = B*B - 4.0*A*C;
  if(discr < 0.0) return -1.0f;
  float temp, sqrt_discrim = sqrtf(discr);
  if (B < 0)
    temp = -0.5f * (B - sqrt_discrim);
  else
    temp = -0.5f * (B + sqrt_discrim);
  const float t0 = temp / A;
  const float t1 = C / temp;

  float t;
  if(t0 <= 0.0f) t = t1;
  else if(t1 <= 0.0f) t = t0;
  else
  {
    t = fminf(t0, t1);
    for(int i=0;i<2;i++)
    {
      for(int k=0;k<3;k++) out[k] = o[k] + t*w[k];
      if(out[0] >= 0.0 && out[0] <= dlen)
        return t;
      t = fmaxf(t0, t1);
    }
    return -1.0f;
  }

  for(int k=0;k<3;k++)
    out[k] = o[k] + t*w[k];
  if(out[0] >= 0.0 && out[0] <= dlen)
    return t;
  else return -1.0f;
}
#endif

static inline float _geo_line_intersect_cone(
    const float *v0,
    const float *v1,
    const float r0,
    const float r1,
    const ray_t *ray,
    float dist,
    hit_t *hit)
{
  float iraylen = 1.0f; // shadow rays aren't normalised:
  if(!hit) iraylen = 1.0f/sqrtf(dotproduct(ray->dir,ray->dir));
  float d[3];
  for(int k=0;k<3;k++) d[k] = v1[k] - v0[k];
  const float d_len = sqrtf(dotproduct(d,d));
  for(int k=0;k<3;k++) d[k] *= 1.0/d_len;
  const float cos_dr = dotproduct(d, ray->dir)*iraylen;
  const float cos_a2 = d_len*d_len/((r1-r0)*(r1-r0) + d_len*d_len);
  // tip of the cone and ray origin wrt tip
  float tip[3], o[3];
  const float tt = -r0*d_len/(r1-r0);
  for(int k=0;k<3;k++) tip[k] = v0[k] + tt*d[k];
  for(int k=0;k<3;k++) o[k] = ray->pos[k] - tip[k];
  const float cos_do = dotproduct(d, o);
  const float cos_ro = dotproduct(ray->dir, o)*iraylen;
  const float cos_oo = dotproduct(o, o);
  const float c2 = cos_dr*cos_dr - cos_a2;
  const float c1 = cos_dr*cos_do - cos_a2*cos_ro;
  const float c0 = cos_do*cos_do - cos_a2*cos_oo;
  // solve the quadratic.  Keep only those X for which Dot(A,X-V) >= 0.
  float tmin = -1.0f;
  if(fabsf(c2) > 0.0)
  {
    const float discr = c1*c1 - c0*c2;
    if(discr < 0.0f) return -1.0f;
    const float root = sqrtf(discr);
    float x[3];
    // go through two solutions, start with closer one:
    for(int i=-1;i<2;i+=2)
    {
      const float t = (-c1 + i*root)/c2;
      if(t > 0.0 && t < dist/iraylen)
      {
        for(int k=0;k<3;k++)
          x[k] = ray->pos[k] + t*ray->dir[k]*iraylen - v0[k];
        const float dt = dotproduct(x, d);
        if(dt >= 0.0f && dt <= d_len)
        {
          if(hit)
          {
            hit->u = dt/d_len;
            float a[3], b[3];
            get_onb(d, a, b);
            hit->v = atan2f(dotproduct(a, x), dotproduct(b, x))/(2.0f*M_PI);
          }
          tmin = dist = t;
        }
      }
    }
  }
  // else some corner cases where the ray is on the cone surface etc.
  return tmin;
}

static inline int geo_line_intersect(
    const prims_t *p,
    const primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const int linestrip = MAX(r0,r1) <= 1e-2f;
  if(linestrip && primid_eq(ray->ignore, pi)) return 0;
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, ray->time, &v0);
  geo_get_vertex_time(p, pi, 1, ray->time, &v1);

  float out[3], len;
  if(fabsf(r1-r0) < 1e-3)
  { // cylinder
    const float t = _geo_line_intersect_cylinder(v0.f, v1.f, r0, ray, out, &len);

    if(t > ray->min_dist && t < hit->dist)
    {
      // init hitpoint
      hit->dist = t;
      hit->prim = pi;
      hit->u = out[0]/len;
      // TODO: consistency check with sample!
      hit->v = atan2f(out[1], out[2])/(2.0f*M_PI);
      return 1;
    }
  }
  else
  { // truncated cone
    const float t = _geo_line_intersect_cone(v0.f, v1.f, r0, r1, ray, hit->dist, hit);
    if((linestrip && t > MAX(ray->min_dist, 1e-3f)) || (!linestrip && t > ray->min_dist))
    {
      hit->dist = t;
      hit->prim = pi;
      return 1;
    }
  }
  return 0;
}


static inline int geo_line_intersect_visible(
    const prims_t *p,
    const primid_t pi,
    const ray_t *ray,
    const float max_dist)
{
  float r0, r1;
  geo_line_get_radii(p, pi, &r0, &r1);
  const int linestrip = MAX(r0,r1) <= 1e-2f;
  if(linestrip && primid_eq(ray->ignore, pi)) return 0;
  float4_t v0, v1;
  geo_get_vertex_time(p, pi, 0, ray->time, &v0);
  geo_get_vertex_time(p, pi, 1, ray->time, &v1);

  if(fabsf(r1-r0) < 1e-3)
  { // cylinder
    float out[3];
    const float t = _geo_line_intersect_cylinder(v0.f, v1.f, r0, ray, out, 0);

    if(((linestrip && t > MAX(ray->min_dist, 1e-3f)) || (!linestrip && t > ray->min_dist)) && t <= max_dist)
      return 1;
  }
  else
  { // truncated cone
    const float t = _geo_line_intersect_cone(v0.f, v1.f, r0, r1, ray, max_dist, 0);
    if(t > ray->min_dist)
      return 1;
  }
  return 0;
}

#endif
