#ifndef _CORONA_SHELL_H
#define _CORONA_SHELL_H

#include "geo/shell_proc.h"

static inline float _geo_shell_extrusion(const prims_t *p, const primid_t pi)
{
  return 2.0f; // XXX read from radius? texture? proc. callback?
}

static inline void geo_shell_get_bounds_shutter_open(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const float *v = geo_get_vertex(p, pi, 0);
  float n[3];
  geo_get_normal_shutter_open(p, pi, 0, n);
  const float e = _geo_shell_extrusion(p, pi);
  *min = fminf(v[dim], v[dim]+e*n[dim]);
  *max = fmaxf(v[dim], v[dim]+e*n[dim]);
  for(int k=1;k<3;k++)
  {
    const float *v = geo_get_vertex(p, pi, k);
    geo_get_normal_shutter_open(p, pi, k, n);
    *min = fminf(v[dim], *min);
    *max = fmaxf(v[dim], *max);
    *min = fminf(v[dim]+e*n[dim], *min);
    *max = fmaxf(v[dim]+e*n[dim], *max);
  }
}

static inline void geo_shell_get_bounds_shutter_close(const prims_t *p, primid_t pi, int dim, float *min, float *max)
{
  const float *v = geo_get_vertex_shutter_close(p, pi, 0);
  float n[3];
  geo_get_normal_shutter_close(p, pi, 0, n);
  const float e = _geo_shell_extrusion(p, pi);
  *min = fminf(v[dim], v[dim]+e*n[dim]);
  *max = fmaxf(v[dim], v[dim]+e*n[dim]);
  for(int k=1;k<3;k++)
  {
    const float *v = geo_get_vertex_shutter_close(p, pi, k);
    geo_get_normal_shutter_close(p, pi, k, n);
    *min = fminf(v[dim], *min);
    *max = fmaxf(v[dim], *max);
    *min = fminf(v[dim]+e*n[dim], *min);
    *max = fmaxf(v[dim]+e*n[dim], *max);
  }
}

static inline float geo_shell_get_area(const prims_t *p, primid_t pi)
{
  // cannot compute area :(
  return 0.0f;
}


static inline void geo_shell_sample(
    const prims_t *p,
    primid_t pi,
    float u,
    float v,
    hit_t *hit,
    const float time)
{
  // cannot sample :(
}

static inline void _geo_shell_st(const prims_t *p, const primid_t pi, const float u, const float v, float *ss, float *tt)
{
  float uv0[3], uv1[3], uv2[3];
  geo_get_uv(p, pi, 0, uv0);
  geo_get_uv(p, pi, 1, uv1);
  geo_get_uv(p, pi, 2, uv2);
  *ss = (1.0f-u-v)*uv0[0] + v*uv1[0] + u*uv2[0];
  *tt = (1.0f-u-v)*uv0[1] + v*uv1[1] + u*uv2[1];
}

// compute normal of triangle
static inline void _geo_shell_get_tri_normal(float v[3][3], float *n)
{
  float edge1[3], edge2[3];
  for(int k=0;k<3;k++)
  {
    edge1[k] = v[1][k] - v[0][k];
    edge2[k] = v[2][k] - v[0][k];
  }
  crossproduct(edge1, edge2, n);
}

// return vertex index 
static inline void _geo_shell_get_vertex(
    const prims_t *p,
    const primid_t pi,
    int ti,            // triangle index, 0-7
    int vi,            // vertex index 0-2
    float time,        // time (motion blur)
    float4_t *v)       // output will be stored here
{
  int i; // triangle index in base mesh [0,3) or extruded version [3,6)
  if(ti == 0)
  { // base triangle
    i = vi;
  }
  else if(ti == 1)
  { // extruded triangle
    i = 3+vi;
  }
  else
  { // sides
    const int side = (ti-2)/2;
    const int si = ti-2 - side*2;
    const int swap = (pi.extra & (1<<side)) >> side;
    int i0 = side+3, i1 = 3+((side+1)%3), i2 = side, i3 = (side+1)%3;
    const int ind[2][2][3] = 
     {{{i0, i2, i1},{i2, i3, i1}}, // regular triangulation
      {{i0, i3, i1},{i2, i3, i0}}}; // swapped
    assert(swap >= 0 && swap < 2);
    assert(si >= 0 && si < 2);
    assert(vi >= 0 && vi < 3);
    i = ind[swap][si][vi];
  }
  assert(i >= 0);
  if(i >= 3) // extruded?
  {
    assert(i < 6);
    float n[3];
    geo_get_vertex_time(p, pi, i-3, time, v);
    geo_get_normal_time(p, pi, i-3, time, n);
    for(int k=0;k<3;k++) v->f[k] += n[k] * _geo_shell_extrusion(p, pi);
  }
  else geo_get_vertex_time(p, pi, i, time, v);
}

// return interpolated triangle
static inline void _geo_shell_get_tri(
    const prims_t *p,
    const primid_t pi,
    const float time, // ray time
    const float w,    // depth
    float v[3][3])    // output triangle vertices
{
  float4_t vb, ve; // bottom and extruded vertices
  for(int k=0;k<3;k++)
  {
    _geo_shell_get_vertex(p, pi, 0, k, time, &vb);
    _geo_shell_get_vertex(p, pi, 1, k, time, &ve);
    for(int j=0;j<3;j++)
      v[k][j] = (1.0f-w)*vb.f[j] + w*ve.f[j];
  }
}

static inline void geo_shell_get_normal_time(
  const prims_t *p,
  primid_t pi,
  hit_t *hit,
  float time)
{
  float n[4];
  float ss, tt;

  _geo_shell_st(p, pi, hit->u, hit->v, &ss, &tt);
  _geo_shell_tex_get_normal(p, pi, ss, tt, hit->w, n);

  // now transform back to world space:
  float tri[3][3];
  float trin[3];
  _geo_shell_get_tri(p, pi, time, hit->w, tri);
  _geo_shell_get_tri_normal(tri, trin);


  // find tangent vectors a and b in world space which correspond to (0,1) and (1,0) in texture (st) space:
  // prims_dmap_t *map = p->maps + t->tex;
  // m * (u,v,1) = (s,t,.) in homogeneous coordinates
  float uv0[2], uv1[2], uv2[2];
  geo_get_uv(p, pi, 0, uv0);
  geo_get_uv(p, pi, 1, uv1);
  geo_get_uv(p, pi, 2, uv2);
  const float m[9] = {uv2[0], uv1[0], uv0[0],
                      uv2[1], uv1[1], uv0[1],
                      1.0f,           1.0f,           1.0f};
  float invm[9];
#define A(y, x) m[(y - 1) * 3 + (x - 1)]
#define B(y, x) invm[(y - 1) * 3 + (x - 1)]
    const float det =
      A(1, 1) * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3)) -
      A(2, 1) * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3)) +
      A(3, 1) * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));
    if(det == 0.0f)
    {
      memcpy(hit->n, trin, sizeof(float)*3);
      memcpy(hit->gn, trin, sizeof(float)*3);
      return; // :(
    }

    const float invDet = 1.f / det;
    B(1, 1) =  invDet * (A(3, 3) * A(2, 2) - A(3, 2) * A(2, 3));
    B(1, 2) = -invDet * (A(3, 3) * A(1, 2) - A(3, 2) * A(1, 3));
    B(1, 3) =  invDet * (A(2, 3) * A(1, 2) - A(2, 2) * A(1, 3));

    B(2, 1) = -invDet * (A(3, 3) * A(2, 1) - A(3, 1) * A(2, 3));
    B(2, 2) =  invDet * (A(3, 3) * A(1, 1) - A(3, 1) * A(1, 3));
    B(2, 3) = -invDet * (A(2, 3) * A(1, 1) - A(2, 1) * A(1, 3));

    B(3, 1) =  invDet * (A(3, 2) * A(2, 1) - A(3, 1) * A(2, 2));
    B(3, 2) = -invDet * (A(3, 2) * A(1, 1) - A(3, 1) * A(1, 2));
    B(3, 3) =  invDet * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2));
#undef A
#undef B
  const float a_st[3] = {ss+1.0f, tt, 1.0f};
  const float b_st[3] = {ss, tt+1.0f, 1.0f};
  float aa[3], ba[3], a[3], b[3];
  memset(ba, 0, sizeof(ba));
  memset(aa, 0, sizeof(aa));
  for(int k=0;k<3;k++)
    for(int l=0;l<3;l++) aa[k] += invm[3*k + l]*a_st[l];
  for(int k=0;k<3;k++)
    for(int l=0;l<3;l++) ba[k] += invm[3*k + l]*b_st[l];
  for(int k=0;k<3;k++) a[k] = aa[0]*tri[2][k] + aa[1]*tri[1][k] + aa[2]*tri[0][k] - hit->x[k];
  for(int k=0;k<3;k++) b[k] = ba[0]*tri[2][k] + ba[1]*tri[1][k] + ba[2]*tri[0][k] - hit->x[k];

  normalise(trin);

  const float depth = _geo_shell_extrusion(p, pi);
  float matrix[9];
  for(int k=0;k<3;k++)
  {
    matrix[k]   = a[k];
    matrix[3+k] = b[k];
    matrix[6+k] = trin[k]/depth; // depth seems to be correctly handled for the unit quad test scene, where st==uv
  }
  // transform with matrix -1 ^T
  hit->n[0] = matrix[0]*n[0] + matrix[3]*n[1] + matrix[6]*n[2];
  hit->n[1] = matrix[1]*n[0] + matrix[4]*n[1] + matrix[7]*n[2];
  hit->n[2] = matrix[2]*n[0] + matrix[5]*n[1] + matrix[8]*n[2];
  normalise(hit->n);
  for(int k=0;k<3;k++) hit->gn[k] = hit->n[k];
}

static inline void geo_shell_init_derivatives(
    const prims_t *p,
    primid_t pi,
    const hit_t *hit,
    const float time,
    vertex_manifold_t *mf)
{
}


static inline float _geo_shell_min_free_path(
    const prims_t *p,
    const primid_t pi,
    const float u,
    const float v,
    const float w,
    const float du,
    const float dv,
    const float dw,
    float *fp) // closest point on free path in uvw space will be returned
{
#if 0 // DEBUG: unit quad has st == uv
  const float d = _prims_tex_min_free_path(p, pi, u, v, w);
  const float len = d/sqrtf(du*du + dv*dv + dw*dw);
  fp[0] = u + len*du;
  fp[1] = v + len*dv;
  fp[2] = w + len*dw;
  return d;
#endif
  // transform (uvw) to texture space, transform (duvw) to texture space (st coordinates).
  // then get min free path and get (uvw + freepath * duvw) in texture space
  // transform that point back and return it.
  float s0, t0, s1, t1;
  _geo_shell_st(p, pi, u, v, &s0, &t0);
  _geo_shell_st(p, pi, u+du, v+dv, &s1, &t1);

  // distance will live on normalized [0,1]^3 coordinates for (s,t,w)
  // TODO: modular backend depending on shapeid!
  const float dist = _geo_shell_tex_min_free_path(p, pi, s0, t0, w);
  // normalise distance in texture space, multiply by free path:
  const float norm = fabsf(dist)/sqrtf((s1-s0)*(s1-s0) + (t1-t0)*(t1-t0) + dw*dw);
  const float ds = norm * (s1-s0);
  const float dt = norm * (t1-t0);
  // const float dw = norm * dw;

  // transform (ds dt dw) back to (uvw) space

  // this matrix is const for every shell and could be cached.
  float uv0[2], uv1[2], uv2[2];
  geo_get_uv(p, pi, 0, uv0);
  geo_get_uv(p, pi, 1, uv1);
  geo_get_uv(p, pi, 2, uv2);
  // m * (du,dv) = (ds,dt)
  float m[4] = {uv2[0] - uv0[0], uv1[0] - uv0[0],
                uv2[1] - uv0[1], uv1[1] - uv0[1]};
  float invm[4];
  const float st[3] = {ds, dt, 1.0f};
  const float det = m[0]*m[3] - m[2]*m[1];
  if(fabsf(det) < 0.001f)
  {
    // degenerate :(
    return 10000.0f;
  }
  invm[0] = m[3]/det;
  invm[3] = m[0]/det;
  invm[1] = -m[1]/det;
  invm[2] = -m[2]/det;

  fp[0] = u + invm[0] * st[0] + invm[1] * st[1];
  fp[1] = v + invm[2] * st[0] + invm[3] * st[1];
  fp[2] = w + norm*dw;
  return dist;
}

// return (unnormalised) normal of convex hull prism.
static inline void _geo_shell_get_normal(
    const prims_t *p,
    const primid_t pi,
    const float time,  // ray time (motion blur)
    const int ti,      // triangle index 0-7
    float *n)          // output: normal
{
  float4_t v0, v1, v2;
  _geo_shell_get_vertex(p, pi, ti, 0, time, &v0);
  _geo_shell_get_vertex(p, pi, ti, 1, time, &v1);
  _geo_shell_get_vertex(p, pi, ti, 2, time, &v2);
  float edge1[3], edge2[3];
  for(int k=0;k<3;k++)
  {
    edge1[k] = v1.f[k] - v0.f[k];
    edge2[k] = v2.f[k] - v0.f[k];
  }
  crossproduct(edge1, edge2, n);
}

// intersect plane of triangle and return uv coords (whether inside or not)
static inline void _geo_shell_intersect_tri(float tri[3][3], const ray_t *ray, float *u, float *v)
{
  float edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  float det, inv_det;

  for(int k=0;k<3;k++)
  {
    edge1[k] = tri[1][k] - tri[0][k];
    edge2[k] = tri[2][k] - tri[0][k];
  }

  crossproduct(ray->dir, edge2, pvec);
  det = dotproduct(edge1, pvec);
  inv_det = 1.0 / det;

  for(int k=0;k<3;k++) tvec[k] = ray->pos[k] - tri[0][k];

  *v = dotproduct(tvec, pvec) * inv_det;
  crossproduct(tvec, edge1, qvec);
  *u = dotproduct(ray->dir, qvec) * inv_det;
}

// get uv coords of point in plane of triangle
static inline void _geo_shell_get_uv(float tri[3][3], float *p, float *u, float *v)
{
  float edge1[3], edge2[3], hit[3];

  for(int k=0;k<3;k++)
  {
    edge1[k] = tri[1][k] - tri[0][k];
    edge2[k] = tri[2][k] - tri[0][k];
    hit[k] = p[k] - tri[0][k];
  }

  float n[3], tmp[3];
  crossproduct(edge1, edge2, n);
  const float len2 = dotproduct(n, n);
  crossproduct(edge1, hit, tmp);
  *u = dotproduct(tmp, n)/len2;
  crossproduct(hit, edge2, tmp);
  *v = dotproduct(tmp, n)/len2;
}

// compute depth coordinate (w) at point of entry p[3].
static inline float _geo_shell_find_height(const prims_t *p, const primid_t pi, const float *x, const float time)
{
  float wM = 1.0f, wm = 0.0f;
  // triangle interpolated for given height w
  float tri[3][3];
  float n[3];
  // TODO: this probably makes sense to optimize for w==1 (enter the shell from above)
  // 12 is from the paper
  float mid;
  for(int i=0;i<12;i++)
  {
    mid = (wM+wm)*.5f;
    _geo_shell_get_tri(p, pi, time, mid, tri);
    _geo_shell_get_tri_normal(tri, n);
    float tmp[3];
    for(int k=0;k<3;k++) tmp[k] = x[k] - tri[0][k];
    if(dotproduct(tmp, n) < 0.0f)
      wM = mid;
    else
      wm = mid;
  }
  return (wM+wm)*.5f;
}

// if the texture-space interpolated point is too far away from the ray (check with dot ray dir), correct:
// construct triangle at height w, intersect ray -> (uvw) back to texture space
static inline void _geo_shell_correct_with_worldspace(
    const prims_t *p,
    const primid_t pi,
    const ray_t *ray,
    float u,
    float v,
    float w,
    float *uu,
    float *vv)
{
  const float tolerance = 0.05f;
  // paper says 4 iterations.
  for(int k=0;k<4;k++)
  {
    float tri[3][3];
    _geo_shell_get_tri(p, pi, ray->time, w, tri);
    float p[3];
    for(int i=0;i<3;i++) p[i] = (1.0f-u-v)*tri[0][i] + v*tri[1][i] + u*tri[2][i] - ray->pos[i];
    const float dot = dotproduct(ray->dir, p);
    if(dot < 1.0f - tolerance)
      // intersect ray with the triangle in world space, try again:
      _geo_shell_intersect_tri(tri, ray, &u, &v);
    else break;
  }
  // return uvs and get out
  *uu = u;
  *vv = v;
}

static inline void geo_shell_intersect(
    const prims_t *p,
    primid_t pi,
    const ray_t *ray,
    hit_t *hit)
{
  // first intersect bounding prism to get entry and exit points.
  float entry_dist = 0;
  float exit_dist  = FLT_MAX;

  // if we are a secondary ray started here, we already know where we're coming from
  // if(hit->prim != pi) // XXX this will be reset to -1 .. :(

  // since we're intersecting a convex hull, it's enough to test all the plane equations (like aabb slabs)
  // and track the max and min of entry and exit points.
  // TODO: sse/avx?
  for(int ti=0;ti<8;ti++)
  {
    float4_t v0;
    float n[3];
    _geo_shell_get_normal(p, pi, ray->time, ti, n);
    _geo_shell_get_vertex(p, pi, ti, 0, ray->time, &v0);

    // flip normal of base triangle to point outside
    const float dot = (ti?1:-1)*dotproduct(ray->dir, n);
    const float dist = (dotproduct(v0.f, n) - dotproduct(ray->pos, n))/dot;
    if(dot < 0.0f) entry_dist = dist > entry_dist ? dist : entry_dist;
    else           exit_dist  = dist < exit_dist  ? dist : exit_dist;
  }
  // early out, could be moved inside the above loop:
  if(exit_dist < 0.0f || entry_dist > hit->dist || exit_dist < entry_dist) return;

  // world space points where we enter and exit:
  float entry[3], exit[3];
  for(int k=0;k<3;k++)
  {
    entry[k] = ray->pos[k] + entry_dist * ray->dir[k];
    exit [k] = ray->pos[k] + exit_dist  * ray->dir[k];
  }

  // find heights w_entry and w_exit
  float w_entry = 0.0f, w_exit = 0.0f;
  w_entry = _geo_shell_find_height(p, pi, entry, ray->time);
  w_exit  = _geo_shell_find_height(p, pi, exit, ray->time);

#if 0
  // DEBUG: show shells
  hit->prim = pi;
  hit->dist = entry_dist;
  hit->u = w_entry;
  hit->v = w_exit;
  hit->w = exit_dist - entry_dist;
  return;
#endif

  float u, v, w, u_entry, v_entry, u_exit, v_exit;
  // get uv for interpolated triangle at w_entry and w_exit of entry and exit points
  float tri[3][3];
  _geo_shell_get_tri(p, pi, ray->time, w_entry, tri);
  _geo_shell_get_uv(tri, entry, &u_entry, &v_entry);
  _geo_shell_get_tri(p, pi, ray->time, w_exit, tri);
  _geo_shell_get_uv(tri, exit, &u_exit, &v_exit);
  u = u_entry; v = v_entry; w = w_entry;

#if 1
  // DEBUG: show uv of entry point
  hit->prim = pi;
  hit->dist = entry_dist;
  hit->u = u_entry;
  hit->v = v_entry;
  hit->w = 0.0f;
  return;
#endif
#if 0
  // DEBUG: show uv of exit point
  hit->prim = pi;
  hit->dist = entry_dist;
  hit->u = u_exit;
  hit->v = v_exit;
  hit->w = 0.0f;
  return;
#endif

  float du, dv, dw;

  du = u_exit - u_entry;
  dv = v_exit - v_entry;
  dw = w_exit - w_entry;
  // going through all the way might be a little much to ask
  du *= 0.1f;
  dv *= 0.1f;
  dw *= 0.1f;

  int count = 0;
  while(1)
  {
    // XXX stupid early out:
    if(count > 10) return;
    // TODO: if dw == 0, this fails!
    // stop with no intersection if exit point is passed
    if((dw > 0.0 && w > w_exit) || (dw < 0.0 && w < w_exit))
        break;

    // step towards exit using the distance map loaded at startup time.
    float fp[3] = {0.0f};
    const float step = _geo_shell_min_free_path(p, pi, u, v, w, du, dv, dw, fp);

    // fprintf(stderr, "[%2d] %f %f %f -> %f %f %f, %f\n", omp_get_thread_num(), u, v, w, du, dv, dw, fp[0]);
    // if an intersection is found (distance map has a 0)
    // check point inside prism (discard intersection if outside)
    if(fabsf(step) < 0.001f)
    {

      // this would be the correct thing to do at ray high curvatures in texture space:
      // if(u >= 0.0f && v >= 0.0f && u+v <= 1.0f)
      {
        float tri[3][3];
        _geo_shell_get_tri(p, pi, ray->time, w, tri);
        for(int i=0;i<3;i++) hit->x[i] = (1.0f-u-v)*tri[0][i] + v*tri[1][i] + u*tri[2][i];
        float x[3];
        for(int i=0;i<3;i++) x[i] = hit->x[i] - ray->pos[i];
        hit->prim = pi;
        hit->dist = dotproduct(ray->dir, x);
#if 0
        hit->u = count/10.;// count to visualize steps
        hit->v = 0;
        hit->w = 0;
#endif
        hit->u = u;
        hit->v = v;
        hit->w = w;
        break;
      }
#if 0
      else
      {
        // make sure we keep stepping through these hits (neighbouring prism, we'll have to ignore it)
        // XXX adjust this to a more sane value! (step = fabsf(step))
        fp[0] = 1./255.0f;
        fp[1] = 1./255.0f;
        fp[2] = 1./255.0f;
      }
#endif
    }

    // try to step along dx by min free path
    const float new_u = fp[0];
    const float new_v = fp[1];
    const float new_w = fp[2];
    float corr_u, corr_v;
    // see if we're totally off track and correct using the world space triangle:
    _geo_shell_correct_with_worldspace(p, pi, ray, new_u, new_v, new_w, &corr_u, &corr_v);
    du = corr_u - u;
    dv = corr_v - v;
    dw = new_w - w;
    // TODO: have a different exit criterion instead, and check the maximum of all stepping directions?
    if(!(fabsf(dw) > 1e-6f)) return;// dw = 0.01f*(w_exit-w_entry); // XXX TODO: this sucks
    // now step that updated direction, using the distance map
    u += du;
    v += dv;
    w += dw;

    count ++;
  }
}

static inline int geo_shell_intersect_visible(
    const prims_t *p,
    primid_t pi,
    const ray_t *ray,
    const float max_dist)
{
  // TODO: optimize for shadow rays!
  hit_t hit;
  hit.dist = max_dist;
  hit.prim = INVALID_PRIMID;
  geo_shell_intersect(p, pi, ray, &hit);
  return !primid_invalid(hit.prim);
}


// tiny helper for shell_init:
static inline void _geo_shell_vert(const prims_t *p, primid_t pi, int i, float4_t *v)
{
  if(i >= 3) // extruded?
  {
    float n[3];
    geo_get_vertex_time(p, pi, i-3, 0.0f, v);
    geo_get_normal_time(p, pi, i-3, 0.0f, n);
    for(int k=0;k<3;k++) v->f[k] += n[k] * _geo_shell_extrusion(p, pi);
  }
  else geo_get_vertex_time(p, pi, i, 0.0f, v);
}

static inline void geo_shell_init(
    prims_t *p,
    uint64_t k)
{
  // fill extra bits with primitive-specific stuff
  // XXX FIXME: motion blur will totally destroy this convex hull property at times!
  p->primid[k].extra = 0;
  for(int side=0;side<3;side++)
  {
    // enumerate sides
    // four vertices on this side:
    int i0 = side+3, i1 = 3+((side+1)%3), i2 = side, i3 = (side+1)%3;
    // try triangulation (i0,i2,i1) and (i2,i3,i1)
    float tmp[3], du[3], dv[3];
    float4_t v0, v1, v2, v3;
    _geo_shell_vert(p, p->primid[k], i0, &v0);
    _geo_shell_vert(p, p->primid[k], i1, &v1);
    _geo_shell_vert(p, p->primid[k], i2, &v2);
    _geo_shell_vert(p, p->primid[k], i3, &v3);
    for(int k=0;k<3;k++) du[k] = v2.f[k] - v0.f[k];
    for(int k=0;k<3;k++) dv[k] = v1.f[k] - v0.f[k];
    crossproduct(du, dv, tmp);
    for(int k=0;k<3;k++) du[k] = v3.f[k] - v0.f[k];
    // swapped order on this side?
    if(dotproduct(tmp, du) > 0.0f)
      p->primid[k].extra |= 1<<side;
  }
  p->primid[k].vcnt = s_prim_shell;
}

#endif
