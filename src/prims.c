#include "prims.h"

#include <assert.h>
#include <string.h>
#include <float.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// slightly bloaty geometry storage backend.

#include "geo.h"
#include "geo/sphere.h"
#include "geo/line.h"
#include "geo/triangle.h"
#include "geo/shell.h"

void prims_get_bounds_shutter_open(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_sphere)
  {
    geo_sphere_get_bounds_shutter_open(p, pi, dim, min, max);
  }
  else if(vcnt == s_prim_line)
  {
    geo_line_get_bounds_shutter_open(p, pi, dim, min, max);
  }
  else if(vcnt == s_prim_shell)
  {
    geo_shell_get_bounds_shutter_open(p, pi, dim, min, max);
  }
  else // if(vcnt == s_prim_tri || vcnt == s_prim_quad)
  {
    geo_tri_get_bounds_shutter_open(p, pi, dim, min, max);
  }
}

void prims_get_bounds_shutter_close(const prims_t *p, primid_t pi, const int dim, float *min, float *max)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_sphere)
  {
    geo_sphere_get_bounds_shutter_close(p, pi, dim, min, max);
  }
  else if(vcnt == s_prim_line)
  {
    geo_line_get_bounds_shutter_close(p, pi, dim, min, max);
  }
  else if(vcnt == s_prim_shell)
  {
    geo_shell_get_bounds_shutter_close(p, pi, dim, min, max);
  }
  else // if(vcnt tri or quad)
  {
    geo_tri_get_bounds_shutter_close(p, pi, dim, min, max);
  }
}

void prims_get_bounds(const prims_t *p, const primid_t pi, int dim, float *min, float *max)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  float m0, M0, m1 = 0, M1 = 0; // init to avoid gcc's stupid maybe-uninitialized logic
  if(vcnt == s_prim_sphere)
  {
    geo_sphere_get_bounds_shutter_open(p, pi, dim, &m0, &M0);
    if(pi.mb)
      geo_sphere_get_bounds_shutter_close(p, pi, dim, &m1, &M1);
  }
  else if(vcnt == s_prim_line)
  {
    geo_line_get_bounds_shutter_open(p, pi, dim, &m0, &M0);
    if(pi.mb)
      geo_line_get_bounds_shutter_close(p, pi, dim, &m1, &M1);
  }
  else if(vcnt == s_prim_shell)
  {
    geo_shell_get_bounds_shutter_open(p, pi, dim, &m0, &M0);
    if(pi.mb)
      geo_shell_get_bounds_shutter_close(p, pi, dim, &m1, &M1);
  }
  else // if(vcnt > 2)
  {
    geo_tri_get_bounds_shutter_open(p, pi, dim, &m0, &M0);
    if(pi.mb)
      geo_tri_get_bounds_shutter_close(p, pi, dim, &m1, &M1);
  }
  if(pi.mb)
  {
    *min = MIN(m0, m1);
    *max = MAX(M0, M1);
  }
  else
  {
    *min = m0;
    *max = M0;
  }
}

void prims_get_aabb(const prims_t *p, primid_t pi, float *aabb)
{
  for(int dim=0;dim<3;dim++)
    prims_get_bounds(p, pi, dim, aabb+dim, aabb+3+dim);
}

void prims_get_aabb_shutter_close(const prims_t *p, primid_t pi, float *aabb)
{
  for(int dim=0;dim<3;dim++)
    prims_get_bounds_shutter_close(p, pi, dim, aabb+dim, aabb+3+dim);
}

void prims_get_shape_aabb(const prims_t *p, uint32_t shapeid, float *aabb)
{
  assert(shapeid < p->num_shapes);
  for(int i=0;i<3;i++) aabb[i] = FLT_MAX;
  for(int i=0;i<3;i++) aabb[3+i] = -FLT_MAX;
  for(uint64_t k=0;k<p->shape[shapeid].num_prims;k++)
  {
    float aabb2[6];
    primid_t pi = p->shape[shapeid].primid[k];
    pi.shapeid = shapeid;
    prims_get_aabb(p, pi, aabb2);
    for(int i=0;i<3;i++)
    {
      aabb[i] = MIN(aabb[i], aabb2[i]);
      aabb[i+3] = MAX(aabb[i+3], aabb2[i+3]);
    }
  }
}

float prims_get_area(const prims_t *p, primid_t pi)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_sphere) return geo_sphere_get_area(p, pi);
  if(vcnt == s_prim_line)   return geo_line_get_area(p, pi);
  if(vcnt == s_prim_shell)  return geo_shell_get_area(p, pi);

  // shutter open area:
  const float *v0 = geo_get_vertex(p, pi, 0);
  const float *v1 = geo_get_vertex(p, pi, 1);
  const float *v2 = geo_get_vertex(p, pi, 2);
  if(vcnt == 3)
    return geo_tri_get_area(v0, v1, v2);
  else if(vcnt == 4)
  {
    const float *v3 = geo_get_vertex(p, pi, 3);
    return geo_tri_get_area(v0, v1, v2) +
      geo_tri_get_area(v0, v2, v3);
  }
  return 0.0f;
}

float prims_get_area_time(const prims_t *p, primid_t pi, const float time)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  // radius isn't motion blurred, area stays:
  if(vcnt == s_prim_sphere) return geo_sphere_get_area(p, pi);
  if(vcnt == s_prim_line)   return geo_line_get_area_time(p, pi, time);
  if(vcnt == s_prim_shell)  return geo_shell_get_area(p, pi); // TODO: time dependent (currently none is implemented)

  float4_t v0, v1, v2, v3;
  geo_get_vertex_time(p, pi, 0, time, &v0);
  geo_get_vertex_time(p, pi, 1, time, &v1);
  geo_get_vertex_time(p, pi, 2, time, &v2);
  if(vcnt == 3)
    return geo_tri_get_area(v0.f, v1.f, v2.f);
  else if(vcnt == 4)
  {
    geo_get_vertex_time(p, pi, 3, time, &v3);
    return geo_tri_get_area(v0.f, v1.f, v2.f) +
      geo_tri_get_area(v0.f, v2.f, v3.f);
  }
  return 0.0f;
}

void prims_retime(const prims_t *p, primid_t pi, hit_t *hit, const float time)
{
  hit->prim = pi;
  float4_t v0, v1, v2, v3;
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_sphere)
  {
    geo_sphere_retime(p, pi, hit, time);
  }
  else if(vcnt == s_prim_line)
  {
    geo_line_retime(p, pi, hit, time);
  }
  else if(vcnt == s_prim_quad)
  {
    geo_get_vertex_time(p, pi, 0, time, &v0);
    geo_get_vertex_time(p, pi, 2, time, &v2);
    if(hit->v >= hit->u)
    {
      geo_get_vertex_time(p, pi, 1, time, &v1);
      geo_tri_retime(v0.f, v1.f, v2.f, hit->u, hit->v-hit->u, hit);
    }
    else
    {
      geo_get_vertex_time(p, pi, 3, time, &v3);
      geo_tri_retime(v0.f, v2.f, v3.f, hit->u-hit->v, hit->v, hit);
    }
  }
  else if(vcnt == s_prim_tri)
  {
    geo_get_vertex_time(p, pi, 0, time, &v0);
    geo_get_vertex_time(p, pi, 1, time, &v1);
    geo_get_vertex_time(p, pi, 2, time, &v2);
    geo_tri_retime(v0.f, v1.f, v2.f, hit->u, hit->v, hit);
  }
  // else // TODO: distance field shell function (how?)
}

void prims_sample(
    const prims_t *p,
    primid_t pi,
    float r0,
    float r1,
    hit_t *hit,
    const float time)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_sphere)
  {
    hit->u = r0;
    hit->v = acosf(r1)/M_PI;
    geo_sphere_retime(p, pi, hit, time);
  }
  else if(vcnt == s_prim_line)
  {
    hit->u = r0; hit->v = r1;
    geo_line_retime(p, pi, hit, time);
  }
  else if(vcnt == s_prim_quad)
  {
    hit->u = r0; hit->v = r1;
    prims_retime(p, pi, hit, time);
  }
  else if(vcnt == s_prim_tri)
  {
    float a = sqrtf(r0);
    float b = (1.0f-r1)*a;
    float c = r1*a;
    a = 1.0f - a;
    hit->u = c;
    hit->v = b;
    prims_retime(p, pi, hit, time);
  }
  // else // TODO: sample distance field shell function (how?)
}

void prims_get_normal_time(const prims_t *p, primid_t pi, hit_t *hit, float time)
{
  const int vcnt = geo_get_vertex_count(p, pi);

  if(vcnt == s_prim_sphere)
  {
    geo_sphere_get_normal_time(p, pi, hit, time);
  }
  else if(vcnt == s_prim_line)
  {
    geo_line_get_normal_time(p, pi, hit, time);
  }
  else if(vcnt == s_prim_shell)
  {
    geo_shell_get_normal_time(p, pi, hit, time);
  }
  else
  {
    // tris and quads:
    float4_t v0, v1, v2, v3;
    float n0[3], n1[3], n2[3], n3[3];
    geo_get_vertex_time(p, pi, 0, time, &v0);
    geo_get_normal_time(p, pi, 0, time, n0);
    geo_get_vertex_time(p, pi, 2, time, &v2);
    geo_get_normal_time(p, pi, 2, time, n2);

    if(vcnt == 3)
    {
      geo_get_vertex_time(p, pi, 1, time, &v1);
      geo_get_normal_time(p, pi, 1, time, n1);
      geo_tri_get_normal(v0.f, v1.f, v2.f, n0, n1, n2, hit->u, hit->v, hit);
    }
    else if(vcnt == 4)
    {
      if(hit->v >= hit->u)
      {
        geo_get_vertex_time(p, pi, 1, time, &v1);
        geo_get_normal_time(p, pi, 1, time, n1);
        geo_tri_get_normal(v0.f, v1.f, v2.f, n0, n1, n2, hit->u, hit->v-hit->u, hit);
      }
      else
      {
        geo_get_vertex_time(p, pi, 3, time, &v3);
        geo_get_normal_time(p, pi, 3, time, n3);
        geo_tri_get_normal(v0.f, v2.f, v3.f, n0, n2, n3, hit->u-hit->v, hit->v, hit);
      }
    }
  }

  // texture coords:
  if(p->shape[pi.shapeid].vtxidx[pi.vi + 0].uv == 0)
  {
    hit->s = hit->u;
    hit->t = hit->v;
  }
  else
  { // have uvs?
    // fill r,s,t coordinates
    hit->r = .0f;
    float uv0[3], uv1[3], uv2[3], uv3[3];
    if(vcnt == s_prim_sphere)
    {
      geo_get_uv(p, pi, 0, uv0);
      hit->s = hit->u + uv0[0];
      hit->t = hit->v + uv0[1];
    }
    else if(vcnt == s_prim_line)
    {
      geo_get_uv(p, pi, 0, uv0);
      geo_get_uv(p, pi, 1, uv1);
      hit->s = uv0[0];
      hit->t = uv0[1];
      hit->r = (1.0f-hit->u) * uv0[2] + hit->u * uv1[2];
    }
    else if(vcnt == s_prim_shell)
    {
      geo_get_uv(p, pi, 0, uv0);
      geo_get_uv(p, pi, 2, uv2);
      geo_get_uv(p, pi, 1, uv1);
      hit->s = (1.0f-hit->u-hit->v)*uv0[0] + hit->v*uv1[0] + hit->u*uv2[0];
      hit->t = (1.0f-hit->u-hit->v)*uv0[1] + hit->v*uv1[1] + hit->u*uv2[1];
    }
    else
    {
      // tris and quads
      geo_get_uv(p, pi, 0, uv0);
      geo_get_uv(p, pi, 2, uv2);
      if(vcnt == 3)
      {
        geo_get_uv(p, pi, 1, uv1);
        hit->s = (1.0f-hit->u-hit->v)*uv0[0] + hit->v*uv1[0] + hit->u*uv2[0];
        hit->t = (1.0f-hit->u-hit->v)*uv0[1] + hit->v*uv1[1] + hit->u*uv2[1];
      }
      if(vcnt == 4)
      {
        if(hit->v >= hit->u)
        {
          geo_get_uv(p, pi, 1, uv1);
          const float u = hit->u, v = hit->v - hit->u;
          hit->s = (1.0f-u-v)*uv0[0] + v*uv1[0] + u*uv2[0];
          hit->t = (1.0f-u-v)*uv0[1] + v*uv1[1] + u*uv2[1];
        }
        else
        {
          geo_get_uv(p, pi, 3, uv3);
          const float u = hit->u-hit->v, v = hit->v;
          hit->s = (1.0f-u-v)*uv0[0] + v*uv2[0] + u*uv3[0];
          hit->t = (1.0f-u-v)*uv0[1] + v*uv2[1] + u*uv3[1];
        }
      }
    }
  }
}

void prims_get_normal(const prims_t *p, primid_t pi, hit_t *hit)
{
  // get normal
  prims_get_normal_time(p, pi, hit, 0.0f);
}

int prims_offset_ray(const hit_t *hit, ray_t *ray)
{
  const float eps = MAX(MAX(.5f, fabsf(hit->x[0])), MAX(fabsf(hit->x[1]), fabsf(hit->x[2])))*1e-4f;

  ray->ignore = hit->prim;
  ray->min_dist = 0.0f;
  // cannot usually do normal offsets, it confuses the hell out of geometric manifolds. we could go two orders of mag
  // closer to zero with eps (1e-5), but directional changes are less forgiving.
  for(int k=0;k<3;k++) ray->pos[k] = hit->x[k] + eps*ray->dir[k];
  // this is sometimes better for spheres and such, but breaks everything
  // that has to do with manifold walks:
  // const float neps = dotproduct(hit->gn, ray->dir) > 0.0f ? 1e-4f : -1e-4f;
  // for(int k=0;k<3;k++) ray->pos[k] = hit->x[k] + eps*ray->dir[k] + neps * hit->gn[k];
  return 0;
}

float prims_get_ray(const hit_t *hit1, const hit_t *hit2, ray_t *ray)
{
  // scale independent epsilon
  const float eps = 1e-4f*MAX(MAX(.5f, fabsf(hit1->x[0])), MAX(fabsf(hit1->x[1]), fabsf(hit1->x[2])));//rt.epsilon;

  ray->ignore = hit1->prim;
  ray->min_dist = 0;

  // find direction
  float dir[3];
  for(int k=0;k<3;k++)
    ray->dir[k] = hit2->x[k] - hit1->x[k];
  const float ilen = 1.0f/sqrtf(dotproduct(ray->dir, ray->dir));
  for(int k=0;k<3;k++) ray->dir[k] *= ilen;

  for(int k=0;k<3;k++)
  {
    if(primid_invalid(hit1->prim))
      ray->pos[k] = hit1->x[k];
    else
      ray->pos[k] = hit1->x[k] + eps * ray->dir[k];
    if(primid_invalid(hit2->prim))
      dir[k] = hit2->x[k] - ray->pos[k];
    else
      dir[k] = hit2->x[k] - eps * ray->dir[k] - ray->pos[k];
  }
  return sqrtf(dotproduct(dir, dir));
}

void prims_init_derivatives(
    const prims_t *p,
    primid_t pi,
    const hit_t *hit,
    const float time,
    vertex_manifold_t *mf)
{
  const int vcnt = geo_get_vertex_count(p, pi);

  if(vcnt == s_prim_sphere)
  {
    // already done when getting dpduv, and only there
    geo_sphere_init_derivatives(p, pi, hit, time, mf);
    return;
  }
  if(vcnt == s_prim_line)
  {
    geo_line_init_derivatives(p, pi, hit, time, mf);
    return;
  }
  else if(vcnt == s_prim_tri || vcnt == s_prim_quad)
  {
    float uv0[2], uv1[2], uv2[2], uv3[2];
    float n0[3], n1[3], n2[3], n3[3];
    float4_t v0, v1, v2, v3;
    if(p->shape[pi.shapeid].vtxidx[pi.vi].uv == -1) goto error;
    geo_get_vertex_time(p, pi, 0, time, &v0);
    geo_get_normal_time(p, pi, 0, time, n0);
    geo_get_vertex_time(p, pi, 2, time, &v2);
    geo_get_normal_time(p, pi, 2, time, n2);
    geo_get_uv(p, pi, 0, uv0);
    geo_get_uv(p, pi, 2, uv2);

    if(vcnt == 3)
    {
      geo_get_vertex_time(p, pi, 1, time, &v1);
      geo_get_normal_time(p, pi, 1, time, n1);
      geo_get_uv(p, pi, 1, uv1);
      geo_tri_dnduv(v0.f, v1.f, v2.f, n0, n1, n2, uv0, uv1, uv2, hit->u, hit->v, hit->a, hit->b, mf->dndu, mf->dndv);
      if(geo_tri_dpduv(v0.f, v1.f, v2.f, uv0, uv1, uv2, hit->u, hit->v, mf->dpdu, mf->dpdv)) goto error;
    }
    else if(vcnt == 4)
    {
      if(hit->v >= hit->u)
      {
        geo_get_vertex_time(p, pi, 1, time, &v1);
        geo_get_normal_time(p, pi, 1, time, n1);
        geo_get_uv(p, pi, 1, uv1);
        geo_tri_dnduv(v0.f, v1.f, v2.f, n0, n1, n2, uv0, uv1, uv2, hit->u, hit->v-hit->u, hit->a, hit->b, mf->dndu, mf->dndv);
        if(geo_tri_dpduv(v0.f, v1.f, v2.f, uv0, uv1, uv2, hit->u, hit->v-hit->u, mf->dpdu, mf->dpdv)) goto error;
      }
      else
      {
        geo_get_vertex_time(p, pi, 3, time, &v3);
        geo_get_normal_time(p, pi, 3, time, n3);
        geo_get_uv(p, pi, 3, uv3);
        geo_tri_dnduv(v0.f, v2.f, v3.f, n0, n2, n3, uv0, uv2, uv3, hit->u-hit->v, hit->v, hit->a, hit->b, mf->dndu, mf->dndv);
        if(geo_tri_dpduv(v0.f, v2.f, v3.f, uv0, uv2, uv3, hit->u-hit->v, hit->v, mf->dpdu, mf->dpdv)) goto error;
      }
    }
  }
  else // TODO: shells don't have derivatives for now :(
  {
error:
    memset(mf->dndu, 0, sizeof(float)*3);
    memset(mf->dndv, 0, sizeof(float)*3);
    get_onb(hit->gn, mf->dpdu, mf->dpdv);
  }
}

/* triangle uvs: w*v0 + v*v1 + u*v2    */
/* intersect a quad:
 *
 * (0,0) v0 ------ v1 (0,1)
 *       | \       |
 *       |   \     |
 *       |     \   |
 *       |       \ |
 * (1,0) v3 ------ v2 (1,1)
 *
 * tri1: v >= u
 * tri (v0 v1 v2) uv => quad uv = (u, v+u)
 * tri (v0 v2 v3) uv => quad uv = (u+v, v)
 * quad uv => tri uv (u', v'-u')
 * quad uv => tri uv (u'-v', v')
 */
void prims_intersect(const prims_t *p, primid_t pi, const ray_t *ray, hit_t *hit)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_tri || vcnt == s_prim_quad)
  {
    // tris and quads
    float4_t v0, v1, v2, v3;
    geo_get_vertex_time(p, pi, 0, ray->time, &v0);
    geo_get_vertex_time(p, pi, 1, ray->time, &v1);
    geo_get_vertex_time(p, pi, 2, ray->time, &v2);
    if(vcnt == 3)
    {
      geo_tri_intersect(v0.f, v1.f, v2.f, pi, ray, hit);
    }
    else if(vcnt == 4)
    {
      if(geo_tri_intersect(v0.f, v1.f, v2.f, pi, ray, hit))
      {
        hit->v += hit->u;
        return;
      }
      geo_get_vertex_time(p, pi, 3, ray->time, &v3);
      if(geo_tri_intersect(v0.f, v2.f, v3.f, pi, ray, hit))
      {
        hit->u += hit->v;
      }
    }
  }
  else if(vcnt == s_prim_sphere)
    geo_sphere_intersect(p, pi, ray, hit);
  else if(vcnt == s_prim_line)
    geo_line_intersect(p, pi, ray, hit);
  else if(vcnt == s_prim_shell)
    geo_shell_intersect(p, pi, ray, hit);
}

int prims_intersect_visible(const prims_t *p, primid_t pi, const ray_t *ray, const float max_dist)
{
  const int vcnt = geo_get_vertex_count(p, pi);
  if(vcnt == s_prim_tri || vcnt == s_prim_quad)
  {
    // tris and quads
    float4_t v0, v1, v2, v3;
    geo_get_vertex_time(p, pi, 0, ray->time, &v0);
    geo_get_vertex_time(p, pi, 1, ray->time, &v1);
    geo_get_vertex_time(p, pi, 2, ray->time, &v2);

    if(vcnt == 3) return geo_tri_intersect_visible(v0.f, v1.f, v2.f, ray, max_dist);
    else if(vcnt == 4)
    {
      if(geo_tri_intersect_visible(v0.f, v1.f, v2.f, ray, max_dist)) return 1;
      geo_get_vertex_time(p, pi, 3, ray->time, &v3);
      if(geo_tri_intersect_visible(v0.f, v2.f, v3.f, ray, max_dist)) return 1;
      return 0;
    }
  }
  else if(vcnt == s_prim_sphere)
    return geo_sphere_intersect_visible(p, pi, ray, max_dist);
  else if(vcnt == s_prim_line)
    return geo_line_intersect_visible(p, pi, ray, max_dist);
  else if(vcnt == s_prim_shell)
    return geo_shell_intersect_visible(p, pi, ray, max_dist);
  return 0;
}

void prims_init(prims_t *p)
{
  memset(p, 0, sizeof(prims_t));
  for(int k=0;k<3;k++) p->ghost_aabb[k] =  FLT_MAX;
  for(int k=3;k<6;k++) p->ghost_aabb[k] = -FLT_MAX;
}

void prims_cleanup(prims_t *p)
{
  for(int k=0;k<p->num_shapes;k++)
  {
    munmap(p->shape[k].data, p->shape[k].data_size);
    if(p->shape[k].fd > 2)
      close(p->shape[k].fd);
  }
  free(p->shape);
  free(p->primid);
  memset(p, 0, sizeof(prims_t));
}

void prims_allocate(prims_t *p, const uint32_t num_shapes)
{
  p->num_shapes = num_shapes;
  p->shape = (prims_shape_t *)common_alloc(16, sizeof(prims_shape_t)*p->num_shapes);
}

void prims_discard_shape(prims_t *p, uint32_t shape)
{
  // only mark as empty to not change shapeids around. prims_allocate_index below
  // will be called right after a sweep over this, so the prims of this shape
  // will just never make it into the global array.
  //
  // i guess we could also deallocate shape resources here.
  p->num_prims -= p->shape[shape].num_prims;
  p->shape[shape].num_prims = 0;
}

// will be called after all shapes have been loaded.
void prims_allocate_index(prims_t *p)
{
  p->primid = (primid_t *)common_alloc(16, sizeof(primid_t)*p->num_prims);
  uint64_t num_loaded_prims = 0;
  for(int shapeid=0;shapeid<p->num_shapes;shapeid++)
  {
    int shell = !strcmp(p->shape[shapeid].tex, "shell");
    for(int64_t k=0;k<p->shape[shapeid].num_prims;k++)
    {
      // copy over to joint array. mmap will forget that we loaded those pages later on, if needed.
      p->primid[num_loaded_prims + k] = p->shape[shapeid].primid[k];
      p->primid[num_loaded_prims + k].shapeid = shapeid;
      if(shell) geo_shell_init(p, num_loaded_prims + k);
    }
    num_loaded_prims += p->shape[shapeid].num_prims;
  }
}

int prims_load_with_flags(
    prims_t *p,
    const char *filename,
    const char *texture,
    const int shader,
    char flags,
    const char *searchpath)
{
  int open_flags = (flags == 'r') ? O_RDONLY : O_RDWR;
  int mmap_flags = (flags == 'r') ? PROT_READ : (PROT_READ|PROT_WRITE);

  const int shapeid = p->num_loaded_shapes;
  p->shape[shapeid].material = shader;
  // store texture on shape
  strncpy(p->shape[shapeid].tex, texture, sizeof(p->shape[shapeid].tex)-1);
  char geoname[1024];
  snprintf(geoname, 1024, "%s.geo", filename);
  p->shape[shapeid].fd = open(geoname, open_flags);
  if((p->shape[shapeid].fd == -1) && searchpath)
  {
    char sn[1024];
    snprintf(sn, 1024, "%s/%s.geo", searchpath, filename); 
    p->shape[shapeid].fd = open(sn, open_flags);
  }
  if(p->shape[shapeid].fd == -1)
  {
    p->num_shapes--;
    fprintf(stderr, "[prims_load] could not load geo `%s'! decreasing shape count to %d.\n", filename, p->num_shapes);
    return 1;
  }
  p->shape[shapeid].data_size = lseek(p->shape[shapeid].fd, 0, SEEK_END);
  lseek(p->shape[shapeid].fd, 0, SEEK_SET);
  common_readahead(p->shape[shapeid].fd, 0, p->shape[shapeid].data_size);
  p->shape[shapeid].data = mmap(0, p->shape[shapeid].data_size, mmap_flags, MAP_SHARED,
                               p->shape[shapeid].fd, 0);
  close(p->shape[shapeid].fd);
  p->shape[shapeid].fd = -1;
  snprintf(p->shape[shapeid].name, 1024, "%s", filename);
  if(p->shape[shapeid].data == (void *)-1)
  {
    perror("[prims_load] mmap");
    p->num_shapes--;
    return 1;
  }

  const prims_header_t *header = (const prims_header_t *)p->shape[shapeid].data;
  if(header->magic != GEO_MAGIC)
  {
    fprintf(stderr, "[prims_load] geo `%s' magic number mismatch!\n", filename);
    p->num_shapes--;
    munmap(p->shape[shapeid].data, p->shape[shapeid].data_size);
    return 1;
  }
  if(header->version != GEO_VERSION)
  {
    fprintf(stderr, "[prims_load] geo `%s' version %d != %d (corona)\n", filename, header->version, GEO_VERSION);
    p->num_shapes--;
    munmap(p->shape[shapeid].data, p->shape[shapeid].data_size);
    return 1;
  }
  p->shape[shapeid].primid = (primid_t *)(header + 1);
  p->shape[shapeid].num_prims = header->num_prims;
  // map our pointers
  p->shape[shapeid].vtxidx = (prims_vtxidx_t *)((uint8_t*)header + header->vtxidx_offset);
  p->shape[shapeid].vtx = (prims_vtx_t *)((uint8_t*)header + header->vertex_offset);

  p->num_prims += header->num_prims;
  p->num_loaded_shapes ++;
  return 0;
}

void prims_extend_ghost_aabb(prims_t *p, const float *aabb)
{
  for(int k=0;k<3;k++) p->ghost_aabb[k] = MIN(p->ghost_aabb[k], aabb[k]);
  for(int k=3;k<6;k++) p->ghost_aabb[k] = MAX(p->ghost_aabb[k], aabb[k]);
}
