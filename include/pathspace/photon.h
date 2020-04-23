#pragma once

#include "lights.h"

typedef struct photon_path_t
{
  uint32_t vi;      // start vertex index
}
photon_path_t;

typedef struct photon_vert_t
{
  float x[3];       // world space position
  uint16_t flags;
  uint16_t mode;
}
photon_vert_t;

typedef struct photon_knn_point_t
{
  uint32_t axis  :  2; // axis
  uint32_t right : 30; // pointer to children
  uint32_t pi;         // index of referenced path in global buffer
  uint32_t vi;         // index of vertex on the path
}
photon_knn_point_t;

typedef struct photon_map_t
{
  photon_knn_point_t *points;

  photon_path_t *paths;          // cached paths
  photon_vert_t *verts;          // vertex pool used by cached paths
  uint32_t max_paths, max_verts; // allocation sizes
  uint64_t num_pathverts;        // actually traced counts, for reasons of sync in one reg
}
photon_map_t;

static photon_vert_t *photon_get_vert(const photon_map_t *s, const int pi, const int vi)
{
  // vertex numbers are real vertex numbers like on path_t
  return s->verts + s->paths[pi].vi + vi;
}

#if 0 // that would be it, but we don't need it
static inline int photon_path_len(const photon_map_t *s, const int pi)
{
  if(pi == s->num_paths - 1)
    return s->num_verts - s->paths[pi].vi;
  return s->paths[pi+1].vi - s->paths[pi].vi;
}
#endif

static inline photon_map_t *photon_init(
    uint32_t num_paths)
{
  photon_map_t *s = calloc(1, sizeof(*s));
  s->max_paths = num_paths;
  // speculate on low average vertex count:
  s->max_verts = s->max_paths * 5;
  s->points = malloc(sizeof(*s->points)*s->max_verts);
  s->paths = malloc(sizeof(*s->paths)*s->max_paths);
  s->verts = malloc(sizeof(*s->verts)*s->max_verts);
  s->num_pathverts = 0;
  fprintf(stderr, "[photon map] allocating %.02f MB for %.02fM photon paths\n", 
      ((sizeof(photon_knn_point_t) + sizeof(photon_vert_t))*s->max_verts + sizeof(photon_path_t)*s->max_paths)/1024.0f/1024.0f, s->max_paths/(1024.0*1024));
  return s;
}

static inline void photon_cleanup(
    photon_map_t *s)
{
  free(s->points);
  free(s->paths);
  free(s->verts);
  free(s);
}

// ================================================================================
//  photons stored in kd tree suitable for knn lookups:
// ================================================================================

static inline void _knn_aabb(const photon_map_t *s, float* aabb, const photon_knn_point_t *points, int si, int ei)
{
  for(int k=0;k<3;k++) aabb[k+3] = aabb[k] = photon_get_vert(s, points[si].pi, points[si].vi)->x[k];
  for (int i = si+1; i < ei; i++)
  {
    assert(points[i].pi >= 0);
    assert(points[i].pi < (s->num_pathverts>>32));
    for(int k=0;k<3;k++) aabb[k]   = MIN(photon_get_vert(s, points[i].pi, points[i].vi)->x[k], aabb[k]);
    for(int k=0;k<3;k++) aabb[k+3] = MAX(photon_get_vert(s, points[i].pi, points[i].vi)->x[k], aabb[k+3]);
  }
}

#define KNN_SWAP(a, b) {photon_knn_point_t temp = points[a]; points[a] = points[b]; points[b] = temp;}

static inline int _knn_split_sort(const photon_map_t *s, photon_knn_point_t *points, int si, int ei, int axis, float d)
{
  int m = 0; 
  while (si < ei)
  {
    while ((si < ei) && (photon_get_vert(s, points[si  ].pi, points[si  ].vi)->x[axis] <  d)) si++;
    while ((si < ei) && (photon_get_vert(s, points[ei-1].pi, points[ei-1].vi)->x[axis] >= d)) ei--;
    if (si < ei-1) KNN_SWAP(si, ei-1);
  }
  m = si;
  return m;
}

static inline void knn_build_rec(const photon_map_t *s, photon_knn_point_t* points, int si, int ei)
{
  if(si == ei) return; // no samples here.
  if(si == ei - 1)
  { // only one sample left
    points[si].right = 0;
    points[si].axis = 3;
    return;
  }

  // find longest axis and split it in the middle
  float bbox[2*3];
  _knn_aabb(s, bbox, points, si, ei);
  float temp[3]; for(int k=0;k<3;k++) temp[k] = bbox[k+3] - bbox[k];
  int axis = 0;
  for(int k=1;k<3;k++) if(temp[k] > temp[axis]) axis = k;
  float d = (bbox[axis] + bbox[3+axis])*0.5f;

  // find the sample closest to the split
  int pIdx = si;
  float dist = fabsf(d - photon_get_vert(s, points[si].pi, points[si].vi)->x[axis]);
  for (int i = si+1; i < ei; i++)
  {
    float npdist = fabsf(d - photon_get_vert(s, points[i].pi, points[i].vi)->x[axis]);
    if (npdist < dist)
    {
      dist = npdist;
      pIdx = i;
    }
  }

  // move the plane to the sample point
  d = photon_get_vert(s, points[pIdx].pi, points[pIdx].vi)->x[axis];
  KNN_SWAP(si, pIdx); // now swap the sample to the beginning of the list
  int m = _knn_split_sort(s, points, si+1, ei, axis, d);
  points[si].right = m;
  points[si].axis = axis;

  if (si+1 < m) knn_build_rec(s, points, si+1, m); // sort the left part
  if (m < ei)   knn_build_rec(s, points, m, ei);   // sort the right part
  if(m >= ei) points[si].right = 0; // the right part is 0
}
#undef KNN_SWAP

// compute squared distance
static inline float _knn_dist(const float *x, const float *y)
{
  float dist = (x[0] - y[0])*(x[0] - y[0]);
  for(int k=1;k<3;k++)
    dist += (x[k] - y[k])*(x[k] - y[k]);
  return dist;
}

// do a nearest neighbour lookup (only one neighbour)
// return hypothetical number of photons counted in given distance.
static inline int photon_knn_find(
    const photon_map_t *s,
    const float *q,             // query position
    const photon_knn_point_t **res_pt,
    float *res_dist,
    const float mdist)          // maximum distance from query point.
{
  const photon_knn_point_t *points = s->points;
  const float mdist2 = mdist*mdist; // actually we're comparing squared distances

  int stack[128];
  int top = 0;
  stack[top] = 0;

  int cnt = 0;
  *res_dist = FLT_MAX;

  while(top >= 0)
  {
    assert(top < 128);
    const photon_knn_point_t *p = points + stack[top];

    const float dist = _knn_dist(photon_get_vert(s, p->pi, p->vi)->x, q);
    if(dist < mdist2)
    {
      cnt ++;
      if(*res_dist > dist)
      { // fill or replace output:
        *res_pt = p;
        *res_dist = dist;
      }
      // mdist2 = dist; // can't cull, need to count photons
    }

    if(p->axis == 3)
    { // check for leaf
      top--;
    }
    else
    {
      int ax = p->axis;
      float d = q[ax] - photon_get_vert(s, p->pi, p->vi)->x[ax];
      if(d < 0)
      { // left side first
        int left = stack[top] + 1;
        if ((p->right != 0) && (left != p->right) && (d*d < mdist2))
          stack[top++] = p->right;
        stack[top] = left;
      }
      else
      { // right side first
        int left = stack[top] + 1;
        if ((left != p->right) && (d*d < mdist2))
          stack[top++] = left;
        if (p->right != 0)
          stack[top] = p->right;
        else
          top--;
      }
    }
  }
  return cnt;
}

static inline int photon_record_path(
    photon_map_t *c,
    const path_t *path,
    const float throughput)
{
  if(throughput > 0 && path->length > 1)
  {
    int plen = path->length;
    int64_t inc = (1ul<<32) | plen;
    uint64_t idx = __sync_fetch_and_add(&c->num_pathverts, inc);
    uint32_t path_idx = idx >> 32;
    if(path_idx >= c->max_paths)
    {
      __sync_fetch_and_add(&c->num_pathverts, -inc);
      return 1; // path buffer full?
    }
    photon_path_t *p = c->paths + path_idx;
    uint32_t vert_idx = idx & 0xfffffffful;
    if(vert_idx + path->length >= c->max_verts)
    {
      __sync_fetch_and_add(&c->num_pathverts, -inc);
      return 1; // vert buffer full?
    }
    p->vi = vert_idx;

    for(int j=0;j<plen;j++)
    { // record point for each vertex
      c->points[vert_idx].vi = j;
      c->points[vert_idx].pi = path_idx;
      for(int k=0;k<3;k++)
        c->verts[vert_idx].x[k] = path->v[j].hit.x[k];
      c->verts[vert_idx].flags = path->v[j].flags;
      c->verts[vert_idx].mode = path->v[j].mode;
      vert_idx++;
    }
    // if((path_idx & 0xff) == 0xff)
    //   fprintf(stderr, "[photon record] %.02f%% paths, avg v/p %f      \r",
    //       100.0*path_idx/(double)c->max_paths, vert_idx/(double)path_idx);
  }
  return 0;
}

// static inline float photon_pdf_path_merge(const photon_map_t *s, path_t *path, const int v);

static inline int photon_path_merge(
    const photon_map_t *s,     // the photon map cache
    path_t *path)              // the path to merge with a photon. vertices past path->length will be overwritten.
{
  // discard specular merge points here
  const int joint = path->length-1;
  if(path->v[joint].material_modes & s_specular) return 1;

  // closest point lookup
  const float rd_rad = 0.01f; // TODO: implement this: restrict to ray differential of path (raydiff_v1 times a couple of manifold matrices)
  const photon_knn_point_t *res_pt = 0;
  float res_dist;
  const int cnt = photon_knn_find(s, path->v[path->length-1].hit.x, &res_pt, &res_dist, rd_rad);
  if(cnt == 0) return 1; // no photons

  const float radius = rd_rad;
  const int pi = res_pt->pi;
  const int vi = res_pt->vi;
  // ray trace the full path
  // const int phlen = photon_path_len(s, pi);
  if(path->length + vi > PATHSPACE_MAX_VERTS) return 1;

  // this throughput will include the bsdf at the joint vertex.
  // we need to evaluate it up front anyways so we get consistent volume stacks
  // since this requires surface reflection/transmission modes to be set.
  double light_path_throughput = 1.0;
  for(int k=vi-1;k>=0;k--)
  {
    const int v = path->length;
    const photon_vert_t *vert = photon_get_vert(s, pi, k);
    for(int i=0;i<3;i++)
      path->v[v].hit.x[i] = vert->x[i];

    // use direction that aligns with path vertices for bsdf at merge point:
    for(int k=0;k<3;k++) path->e[v].omega[k] = path->v[v].hit.x[k] - path->v[v-1].hit.x[k];
    path->e[v].dist = sqrtf(dotproduct(path->e[v].omega, path->e[v].omega));
    for(int k=0;k<3;k++) path->e[v].omega[k] /= path->e[v].dist;

    // init the material mode flags for meaningful volume stacks:
    const float bsdf = shader_brdf(path, v-1);
    if(bsdf <= 0.0) return 1;

    path->v[v].flags = vert->flags;
    path->v[v].mode = vert->mode;
    if(path_project(path, v, s_propagate_reconstruct)) return 1;
    path->v[v].tech = s_tech_extend_adj;
    if((path->v[v].flags & s_environment) && k)
      return 1; // premature end of path

    // we buy that this is a good vertex now:
    path->length++;

    if(k==vi-1)
    { // support specular paths:
      // set omega[v] to exactly the direction from the cached photon path.
      // this is only needed the first time around, the rest of the vertex
      // positions and directions will be in line with what is stored in the
      // photon map (this is asserted in the end of path_project).
      for(int k=0;k<3;k++) path->e[v].omega[k] = vert->x[k] - path->v[v-1].hit.x[k];
      path->e[v].dist = sqrtf(dotproduct(path->e[v].omega, path->e[v].omega));
      for(int k=0;k<3;k++) path->e[v].omega[k] /= path->e[v].dist;
      // TODO: if envmap, need to re-evaluate shading.em for the correct direction
    }

    light_path_throughput *= bsdf * path_G(path, v) * shader_vol_transmittance(path, v);
    // XXX assumes we start at the sensor!
    if(!k) light_path_throughput *= lights_eval_vertex(path, v);
  }

  // check end point is light or sensor:
  if(!(path->v[path->length-1].mode & (s_emit|s_sensor)))
    return 1; // missed light source or wrong wavelength or something.

  path->v[joint].tech = s_tech_merged;
  for(int k=joint+1;k<path->length;k++)
  {
    path->v[k].pdf = path_pdf_extend_adjoint(path, k);
    light_path_throughput /= path->v[k].pdf;
  }

  const double joint_pdf_adj = path_pdf_extend_adjoint(path, joint);
  light_path_throughput /= joint_pdf_adj;

  // the last two vertices are fine, that points to broken volume nesting,
  // propagated all the way to the end.
  for(int k=0;k<path->length-2;k++)
    if(!(path->v[k].pdf > 0.0f)) return 1;
  if(path->v[path->length-2].pdf <= 0.0f) return 1;

  const float domain = (path->v[joint].mode & s_volume) ?
      4.0/3.0 * radius*radius*radius * M_PI :
      radius*radius * M_PI;
  path->throughput = path->v[joint].throughput * light_path_throughput * cnt/(s->max_paths * domain);
  // path->v[joint].pdf = path->v[joint].pdf * joint_pdf_adj * s->max_paths * domain/cnt;
  path->v[joint].pdf = path->v[joint].pdf * joint_pdf_adj * s->max_paths * domain;

  if(!(path->v[joint].pdf > 0.0)) return 1; // numerical drop to 0
  return 0;
}

// returns vertex area pdf of the given vertex v on the path, as if it had
// been constructed via photon_path_merge().
static inline float photon_pdf_path_merge(
    const photon_map_t *s, // photon map
    path_t *path,          // complete path
    const int v)           // merged vertex
{
  if(path->v[v].mode & s_specular) return 0.0f;
  // return fake pdf that will amount to some (biased) estimator if computing measurement/pdf.
  float pdf = path_pdf_extend(path, v);
  pdf *= path_pdf_extend_adjoint(path, v);
  // unfortunately need closest point lookup to obtain correct radius
  const float rd_rad = 0.01f; // TODO: implement this: restrict to ray differential of path (raydiff_v1 times a couple of manifold matrices)
#if 0
  const photon_knn_point_t *res_pt = 0;
  float res_dist;
  const int cnt = photon_knn_find(s, path->v[v].hit.x, &res_pt, &res_dist, rd_rad);
  if(cnt == 0) return 0.0f; // no photons
#endif

  const float radius = rd_rad;
  const float domain = (path->v[v].mode & s_volume) ?
      4.0/3.0 * radius*radius*radius * M_PI :
      radius*radius * M_PI;
  // pdf *= s->max_paths * domain/cnt;
  pdf *= s->max_paths * domain;
  return pdf;
}

static inline void photon_clear(
    photon_map_t *s)
{
  s->num_pathverts = 0ul;
}

static inline void photon_build(
    photon_map_t *s)
{
  knn_build_rec(s, s->points, 0, s->num_pathverts & 0xfffffffful);
}
