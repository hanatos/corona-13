#ifndef CORONA_GEO_SPHERE_H
#define CORONA_GEO_SPHERE_H

#include "corona_common.h"
#include "sampler_common.h"
#include "pathspace.h"


static inline float geo_sphere_get_radius(const prims_t *p, primid_t pi)
{
  prims_shape_t *s = p->shape + pi.shapeid;
  return tofloat(s->vtx[(pi.mb+1)*s->vtxidx[pi.vi].v].n);
}


static inline void geo_sphere_get_bounds_shutter_open(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const float *v = geo_get_vertex(p, pi, 0);
  const float radius = geo_sphere_get_radius(p, pi);
  *min = v[dim] - radius;
  *max = v[dim] + radius;
}

static inline void geo_sphere_get_bounds_shutter_close(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const float *v = geo_get_vertex_shutter_close(p, pi, 0);
  const float radius = geo_sphere_get_radius(p, pi);
  *min = v[dim] - radius;
  *max = v[dim] + radius;
}

static inline float geo_sphere_get_area(const prims_t *p, primid_t pi)
{
  const float r = geo_sphere_get_radius(p, pi);
  return 4.0f*M_PI*r*r;
}

static inline void geo_sphere_retime(
    const prims_t *p,
    primid_t pi,
    hit_t *hit,
    const float time)
{
  const float r = geo_sphere_get_radius(p, pi);
  sample_sphere(hit->x, hit->x+1, hit->x+2, -(cosf(hit->v*M_PI)-1.f)/2.f, hit->u);
  float4_t c;
  geo_get_vertex_time(p, pi, 0, time, &c);
  for(int k=0;k<3;k++) hit->x[k] = c.f[k] + r*hit->x[k];
}

static inline void geo_sphere_get_normal_time(
    const prims_t *p,
    primid_t pi,
    hit_t *hit,
    const float time)
{
  float4_t c;
  geo_get_vertex_time(p, pi, 0, time, &c);
  for(int k=0;k<3;k++) hit->gn[k] = hit->x[k] - c.f[k];
  normalise(hit->gn);
  memcpy(hit->n, hit->gn, sizeof(float)*3);
}

static inline void sphere_init_derivatives(
    const float *c, // center point
    const float r,  // radius
    const float *p, // position
    vertex_manifold_t *mf)
{
  const float x[3] = {p[0]-c[0], p[1]-c[1], p[2]-c[2]};
  mf->dpdu[0] = -x[1] * 2*M_PI;
  mf->dpdu[1] =  x[0] * 2*M_PI;
  mf->dpdu[2] =  0.0f;
  const float zrad = x[0]*x[0] +  x[1]*x[1];
  const float cos_theta = x[2]/r;
  const float sin_theta = sqrtf(fmaxf(0.0f, 1.0f-cos_theta*cos_theta));
  float cos_phi, sin_phi;
  if(zrad > 0.0f)
  {
    cos_phi = x[0] / zrad;
    sin_phi = x[1] / zrad;
  }
  else
  {
    cos_phi = 0.0f;
    sin_phi = 1.0f;
  }
  mf->dpdv[0] = x[2] * cos_phi * M_PI;
  mf->dpdv[1] = x[2] * sin_phi * M_PI;
  mf->dpdv[2] = -sin_theta * r * M_PI;

  for(int k=0;k<3;k++)
  {
    mf->dndu[k] = mf->dpdu[k] / r;
    mf->dndv[k] = mf->dpdv[k] / r;
  }
}

static inline void geo_sphere_init_derivatives(
    const prims_t *p,
    primid_t pi,
    const hit_t *hit,
    const float time,
    vertex_manifold_t *mf)
{
  float4_t c;
  geo_get_vertex_time(p, pi, 0, time, &c);
  const float r = geo_sphere_get_radius(p, pi);
  sphere_init_derivatives(c.f, r, hit->x, mf);
}

static inline float _geo_sphere_intersect(
    const float *center,
    const float radius,
    const ray_t *ray)
{
  const float a = dotproduct(ray->dir, ray->dir);
  float o[3] = {ray->pos[0]-center[0],ray->pos[1]-center[1],ray->pos[2]-center[2]};
  const float b = 2.0f*dotproduct(o, ray->dir);
  const float c = dotproduct(o, o) - radius*radius;
  if (a == 0)
  {
    if (b != 0)
      return -c / b;
    return -FLT_MAX;
  }

  float discrim = b*b - 4.0f*a*c;
  if (discrim < 0) return -FLT_MAX;

  float temp, sqrt_discrim = sqrtf(discrim);
  if (b < 0)
    temp = -0.5f * (b - sqrt_discrim);
  else
    temp = -0.5f * (b + sqrt_discrim);

  const float x0 = temp / a;
  const float x1 = c / temp;

  // return minimum if not smaller than zero
  if(x0 <= 0.0f) return x1;
  else if(x1 <= 0.0f) return x0;
  else return fminf(x0, x1);
}

static inline int geo_sphere_intersect(
    const prims_t *p,
    primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  float4_t center;
  geo_get_vertex_time(p, pi, 0, ray->time, &center);
  const float radius = geo_sphere_get_radius(p, pi);
  const float t = _geo_sphere_intersect(center.f, radius, ray);
  if(t > ray->min_dist && t < hit->dist)
  {
    hit->dist = t;
    hit->prim = pi;
    for(int k=0;k<3;k++) hit->x[k] = ray->pos[k] + t*ray->dir[k]; // cheaper than trig in get normal
    hit->u = atan2f((hit->x[1]-center.f[1])/radius, (hit->x[0]-center.f[0])/radius)/(2.0f*M_PI);
    hit->v = acosf(CLAMP((hit->x[2]-center.f[2])/radius, -1.0f, 1.0f))/M_PI;
    return 1;
  }
  return 0;
}

static inline int geo_sphere_intersect_visible(
    const prims_t *p,
    primid_t pi,
    const ray_t *ray,
    const float max_dist)
{
  float4_t center;
  geo_get_vertex_time(p, pi, 0, ray->time, &center);
  const float radius = geo_sphere_get_radius(p, pi);
  const float t = _geo_sphere_intersect(center.f, radius, ray);
  if(t > 0.0f && t <= max_dist) return 1;
  return 0;
}

#endif
