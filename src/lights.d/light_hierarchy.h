#pragma once

#include "corona_common.h"

// light hierarchy initially written by florian simon
// and then updated by addis dittebrandt
// stolen from the quake2 code base and adjusted to work here.

typedef struct lh_cone_t
{
  float axis[3]; // mean direction of lights in node
  float th_o;    // max normal deviaton
  float th_e;    // max emission angle
}
lh_cone_t;

typedef struct lh_prim_t
{
  primid_t primid;
  float c[3];
  float energy;
  float aabb[6];
  lh_cone_t cone;
}
lh_prim_t;

typedef struct lh_node_t
{
  uint64_t child[4];     // child pointer or leaf info
  float    energy[4];
  float    aabb[4][6];
  uint32_t cone_axis[4];
  uint8_t  cone_th_o[4];
  uint8_t  cone_th_e[4];
  // 168 bytes, worth 21 uint64_t
}
lh_node_t;

typedef struct light_hierarchy_t
{
  float aabb[6];
  uint64_t max_num_nodes;
  uint64_t num_nodes;
  lh_node_t *nodes;
  uint64_t max_num_prims;
  uint64_t num_prims;
  lh_prim_t *prims;
}
light_hierarchy_t;

// TODO: consider edf/cosine power
static inline lh_cone_t
lh_triangle_to_cone(float* t0, float* t1, float* t2)
{
  lh_cone_t r;
  float u[3], v[3];
  for(int k=0;k<3;k++) u[k] = t1[k] - t0[k];
  for(int k=0;k<3;k++) v[k] = t2[k] - t0[k];
  crossproduct(u, v, r.axis);
  normalise(r.axis);
  r.th_o = 0.0f;
  r.th_e = M_PI/2.0f;
  return r;
}

static inline void
lh_cleanup(light_hierarchy_t* lh)
{
  free(lh->nodes);
  free(lh->prims);
  lh->nodes = 0;
  lh->prims = 0;
}

static inline void
lh_init(light_hierarchy_t *lh)
{
  memset(lh, 0, sizeof(*lh));
  for(int k=0;k<3;k++)
  {
    lh->aabb[k]   =  FLT_MAX;
    lh->aabb[3+k] = -FLT_MAX;
  }
}

static inline lh_cone_t
lh_cone_union(lh_cone_t a, lh_cone_t b)
{
  // don't expand by empty cones
  if(a.axis[0] == 0.0f && a.axis[1] == 0.0f && a.axis[2] == 0.0f)
    return b;
  if(b.axis[0] == 0.0f && b.axis[1] == 0.0f && b.axis[2] == 0.0f)
    return a;

  lh_cone_t r;

  if (b.th_o > a.th_o)
  {
    lh_cone_t t = b;
    b = a;
    a = t;
  }

  float th_d = acosf(CLAMP(dotproduct(a.axis, b.axis), -1.0f, 1.0f));
  r.th_o = (a.th_o + th_d + b.th_o) / 2.0f;
  r.th_e = fmaxf(a.th_e, b.th_e);

  if (fminf(th_d + b.th_o, M_PI) <= a.th_o)
  {
    // a already covers b
    r.th_o = a.th_o;
    for(int k=0;k<3;k++) a.axis[k] = r.axis[k];
  }
  else if (M_PI <= r.th_o)
  {
    r.th_o = M_PI;
    for(int k=0;k<3;k++) a.axis[k] = r.axis[k];
  }
  else
  {
    float th_r = r.th_o - a.th_o;
    float axis_r[3];
    crossproduct(a.axis, b.axis, axis_r);
    normalise(axis_r);
    float *dir = axis_r;
    float *point = a.axis;
    float matrix[3][3];
    float angle, s, c, one_c, xx, yy, zz, xy, yz, zx, xs, ys, zs;

    angle = th_r;
    s = sinf(angle);
    c = cosf(angle);
    one_c = 1.0f - c;

    xx = dir[0] * dir[0];
    yy = dir[1] * dir[1];
    zz = dir[2] * dir[2];
    xy = dir[0] * dir[1];
    yz = dir[1] * dir[2];
    zx = dir[2] * dir[0];
    xs = dir[0] * s;
    ys = dir[1] * s;
    zs = dir[2] * s;

    matrix[0][0] = (one_c * xx) + c;
    matrix[0][1] = (one_c * xy) - zs;
    matrix[0][2] = (one_c * zx) + ys;

    matrix[1][0] = (one_c * xy) + zs;
    matrix[1][1] = (one_c * yy) + c;
    matrix[1][2] = (one_c * yz) - xs;

    matrix[2][0] = (one_c * zx) - ys;
    matrix[2][1] = (one_c * yz) + xs;
    matrix[2][2] = (one_c * zz) + c;

    r.axis[0] = dotproduct(matrix[0], point);
    r.axis[1] = dotproduct(matrix[1], point);
    r.axis[2] = dotproduct(matrix[2], point);
    normalise(r.axis);
  }

  return r;
}

// used during tree construction, bin information
// finding a best split
typedef struct lh_bin_t
{
  int num_prims;
  float c_aabb[6];
  float energy;
  float aabb[6];
  lh_cone_t cone;
}
lh_bin_t;

static inline void
lh_len(float aabb[6], float lengths[3])
{
  for (int k = 0; k < 3; k++)
    lengths[k] = aabb[k + 3] - aabb[k];
}

static inline float
lh_sur_m(float lengths[3])
{
  float area = 0.0f;
  for (int k = 0; k < 3; k++)
    area += 2.0f * lengths[k] * lengths[(k + 1) % 3];
  return area;
}

static inline float
lh_ori_m(lh_cone_t c)
{
  float th_w = fminf(c.th_o + c.th_e, M_PI);
  return 2.0f * M_PI * (1.0f - cosf(c.th_o)) + 0.5f * M_PI * (
      2.0f * th_w * sinf(c.th_o)
      - cosf(c.th_o - 2.0f * th_w)
      - 2.0f * c.th_o * sinf(c.th_o)
      + cosf(c.th_o)
      );
}

static inline float
lh_reg_m(float lengths[3], int dim)
{
  float max_length = fmaxf(fmaxf(lengths[0], lengths[1]), lengths[2]);
  return max_length / lengths[dim];
}

static inline float
lh_cost_measure(
    int dim,
    lh_bin_t p_bin,
    int num_children,
    lh_bin_t bins[])
{
  float p_lengths[3];
  lh_len(p_bin.aabb, p_lengths);

  float cost = 0.0f;
  for (int k = 0; k < num_children; k++)
  {
    if (bins[k].num_prims == 0) continue;
    float lengths[3];
    lh_len(bins[k].aabb, lengths);
    cost += bins[k].energy * lh_sur_m(lengths) * lh_ori_m(bins[k].cone);
  }

  return lh_reg_m(p_lengths, dim) * cost / (lh_sur_m(p_lengths) * lh_ori_m(p_bin.cone));
}

static inline void
lh_init_aabb(float aabb[6])
{
  for (int k = 0; k < 3; k++)
  {
    aabb[k] = FLT_MAX;
    aabb[k + 3] = -FLT_MAX;
  }
}

static inline void
lh_enlarge_aabb_point(float aabb[6], float p[3])
{
  for (int k = 0; k < 3; k++)
  {
    aabb[k] = fminf(aabb[k], p[k]);
    aabb[k+3] = fmaxf(aabb[k+3], p[k]);
  }
}

static inline void
lh_enlarge_aabb_points(float aabb[6], float *p[], int num_p)
{
  for (int j = 0; j < num_p; j++)
    lh_enlarge_aabb_point(aabb, p[j]);
}

static inline void
lh_enlarge_aabb_aabb(float aabb[6], float a_aabb[6])
{
  for (int k = 0; k < 3; k++)
  {
    aabb[k] = fminf(aabb[k], a_aabb[k]);
    aabb[k+3] = fmaxf(aabb[k+3], a_aabb[k+3]);
  }
}

// returns split index between left and right
static inline uint64_t
lh_split(
    lh_prim_t  *prims,
    uint64_t    left,
    uint64_t    right,
    const float aabb[6],
    float       c_aabb[2][6])
{
  float k_0[3], k_1[3];
  int skip[3] = {0};
  for(int d=0;d<3;d++)
  {
    k_0[d] = aabb[d];
    float diff = aabb[d+3] - aabb[d];
    if(diff == 0) skip[d] = 0;
    k_1[d] = (float)num_bins  / (aabb[d+3] - aabb[d]);
  }

  assert(skip[0] || skip[1] || skip[2]);

  const int num_bins = 8;
  lh_bin_t bins[3][num_bins];
  uint64_t num_prims = right - left;
  // cap bin count to prim count
  num_bins = min(num_bins, num_prims);

  // initialize bins
  for(int d=0;d<3;d++) for (int b = 0; b < num_bins; b++)
  {
    lh_bin_t *bin = &bins[d][b];
    bin->num_prims = 0;
    lh_init_aabb(bin->c_aabb);
    bin->energy = 0.0f;
    lh_init_aabb(bin->aabb);
  }

  // populate bins, walk through prims only once for all three axes
  for (uint64_t i = left; i < right; i++)
  {
    for(int d=0;d<3;d++)
    {
      if(skip[d]) continue;
      lh_prim_t *prim = &prims[i];
      int bin_i = min(k_1[d] * (prim->c[d] - k_0), num_bins - 1);
      lh_bin_t *bin = &bins[d][bin_i];

      bin->num_prims++;
      lh_enlarge_aabb_point(bin->c_aabb, prim->c);
      bin->energy += prim->energy;
      lh_enlarge_aabb_aabb(bin->aabb, prim->aabb);
      if (bin->num_prims == 1)
        bin->cone = prim->cone;
      else
        bin->cone = lh_cone_union(bin->cone, prim->cone);
    }
  }

  // initialize accumulated bins
  lh_bin_t a_bins[2][3][num_bins];
  for(int d=0;d<3;d++) for (int i = 0; i < 2; i++)
    memcpy(a_bins[i][d], bins[d], num_bins * sizeof(lh_bin_t));

  // accumulate bins from left and right
  for(int d=0;d<3;d++)
  {
    if(skip[d]) continue;
    for (int s = 0; s < 2; s++)
    {
      for (int b = 1; b < num_bins; b++)
      {
        int b_i = s == 0 ? b : num_bins - b - 1;
        int prev_b_i = s == 0 ? b_i - 1 : b_i + 1;

        lh_bin_t *prev_bin = &a_bins[s][d][prev_b_i];
        lh_bin_t *bin = &a_bins[s][d][b_i];

        if (bin->num_prims == 0)
          memcpy(bin, prev_bin, sizeof(lh_bin_t));
        else
        {
          bin->num_prims += prev_bin->num_prims;
          lh_enlarge_aabb_aabb(bin->c_aabb, prev_bin->c_aabb);
          bin->energy += prev_bin->energy;
          lh_enlarge_aabb_aabb(bin->aabb, prev_bin->aabb);
          bin->cone = lh_cone_union(bin->cone, prev_bin->cone);
        }
      }
    }
  }

  // find best split candidate
  float m_cost = FLT_MAX;
  int m_s = 0, m_d = 0;

  for(int d=0;d<3;d++)
  {
    if(skip[d]) continue;
    lh_bin_t p_bin = a_bins[1][d][0];
    for (int s = 0; s < num_bins - 1; s++)
    {
      lh_bin_t a_bin[] = {a_bins[0][d][s], a_bins[1][d][s + 1]};
      assert(a_bin[0].num_prims != 0 && a_bin[1].num_prims != 0);
      float cost = lh_cost_measure(d, p_bin, 2, a_bin);

      if (cost < m_cost)
      {
        m_cost = cost;
        m_s = s;
        m_d = d;
      }
    }
  }

  assert(m_cost != FLT_MAX);

  int d = m_d;
  int s = m_s;

  lh_bin_t *a_bin[] = {&a_bins[d][0][s], &a_bins[d][1][s + 1]};

  // rearrange primitives among subtrees
  uint64_t left_start = left;
  uint64_t left_end = left + a_bin[0]->num_prims;
  if(left_end == right)
  { // stupid mid-split
    assert(0); // TODO: compute c_aabb[]
    return (left+right)/2;
  }
  uint64_t right_start = left_end;
  uint64_t right_end = right;
  uint64_t l = left_start;
  uint64_t r = right_start;
  uint64_t term = 0;
  while (1)
  {
    for (; l < left_end; l++)
    {
      lh_prim_t *prim = &prims[l];
      int bin_i = k_1[d] * (prim->c[d] - k_0[d]);
      if (s < bin_i)
        break;
    }

    if (l == left_end)
      term = 1;

    for (; r < right_end; r++)
    {
      lh_prim_t *prim = &prims[r];
      int bin_i = k_1[d] * (prim->c[d] - k_0[d]);
      if (bin_i <= s)
        break;
    }

    assert(term ? r == right_end : r < right_end);

    if (term)
      break;

    lh_prim_t t = prims[l];
    prims[l] = prims[r];
    prims[r] = t;
  }

  for (int s = 0; s < 2; s++)
    memcpy(c_aabb[s], a_bin[s]->c_aabb, sizeof(c_aabb[s]));
  return left_end;
}

static inline uint8_t
f2bf(float frac)
{
  return (uint8_t) floor(256.0f * frac + .5f);
}

static inline float
v2f(float value, float min, float max)
{
  if (min == max) return 0.0f;
  else return (value - min) / (max - min);
}

static inline float
b2v(uint8_t byte, float min, float max)
{
  float f = byte / 256.0f;
  return (1 - f) * min + f * max;
  return min + (max - min) * (byte / 256.0f);
}

// this is the actual work in building the tree:
static inline void
lh_build_binned_rec(
    lh_node_t     *parent,
    const int      childid,
    const uint64_t max_num_nodes,
    uint64_t      *num_nodes,
    lh_node_t     *nodes,
    const uint64_t left,
    const uint64_t right,
    lh_prim_t     *prims,
    const float    aabb[6],
    lh_cone_t     *ccone)
{
  if(parent && (right-left <= 1))
  { // only one primitive left, build a leaf:
    uint64_t node_i = 1ul<<63;  // leaf marker
    node_i |= (right-left)<<55; // primitive count, currently one or empty leaf
    node_i |= left;             // primitive start index
    parent->child[childid] = node_i;
    if(right - left)
      *ccone = prims[left].cone;
    parent->cone_axis[childid] = geo_encode_normal(ccone->axis);
    parent->cone_th_o[childid] = f2bf(v2f(ccone->th_o, 0, M_PI));
    parent->cone_th_e[childid] = f2bf(v2f(ccone->th_e, 0, M_PI));
    parent->energy[childid] = (right-left) ? prims[left].energy : 0.0f;
    return;
  }

  // if we're not a leaf, we need to allocate a node:
  lh_node_t *node = nodes[(*num_nodes)++];
  assert(*num_nodes <= max_num_nodes); // we always have enough memory, right?

  // split three times for 4 children:
  float aabb0[2][6], aabb1[2][6], aabb2[2][6];
  uint64_t end0 = lh_split(prims, left, right, aabb, aabb0);
  uint64_t end1 = lh_split(prims, left, end0,  aabb, aabb1);
  uint64_t end2 = lh_split(prims, end0, right, aabb, aabb2);

  // remember child boxes in our node
  memcpy(node->aabb[0], aabb1[0], sizeof(float)*6);
  memcpy(node->aabb[1], aabb1[1], sizeof(float)*6);
  memcpy(node->aabb[2], aabb2[0], sizeof(float)*6);
  memcpy(node->aabb[3], aabb2[1], sizeof(float)*6);

  lh_cone_t cone[4] = {{{0}}};

  // recurse or make leaf?
  lh_build_binned_rec(node, 0, max_num_nodes, num_nodes, nodes, left, end1,  prims, aabb1[0], cone+0);
  lh_build_binned_rec(node, 1, max_num_nodes, num_nodes, nodes, end1, end0,  prims, aabb1[1], cone+1);
  lh_build_binned_rec(node, 2, max_num_nodes, num_nodes, nodes, end0, end2,  prims, aabb2[0], cone+2);
  lh_build_binned_rec(node, 3, max_num_nodes, num_nodes, nodes, end2, right, prims, aabb2[1], cone+3);

  if(parent)
  { // node aggregated from children:
    parent->energy[childid] = node->energy[0];
    *ccone = cone[0];
    for(int k=1;k<4;k++)
    {
      parent->energy[childid] += node->energy[k];
      *ccone = lh_cone_union(*ccone, cone[k]);
    }
    parent->child[childid] = node - nodes;
    parent->cone_axis[childid] = geo_encode_normal(ccone.axis);
    parent->cone_th_o[childid] = f2bf(v2f(ccone.th_o, 0, M_PI));
    parent->cone_th_e[childid] = f2bf(v2f(ccone.th_e, 0, M_PI));
    // parent takes care of its own bounding boxes
  }
}


static inline void
lh_build_binned(light_hierarchy_t *lh)
{
  // primitives with aabb and bounding cones have been inited from the outside
  lh->max_num_nodes = lh->num_prims*2;
  lh->nodes = malloc(sizeof(lh_node_t)*lh->max_num_nodes);
  lh_build_binned_rec(0, 0, lh->max_num_nodes, 0, lh->nodes, 0, lh->num_prims, lh->prims, lh->aabb[0], 0);
}





// ==============================================================
//  run time things: sampling and pdf eval:
// ==============================================================


// computes the spectral pdf of next event estimation for p->v[v]
mf_t lights_pdf_next_event(const path_t *p, int v)
{
  // XXX FIXME TODO
  // traverse light hierarchy and ask for pdf
  return mf_set1(1.0f);
}



static inline float
lh_projected_tri_area(
    lh_prim_t *lhprim,
    float *p, float *n)
{
  float prims_get_area(const prims_t *p, primid_t pi);
  void prims_get_normal(const prims_t *p, primid_t pi, hit_t *hit);
  float a[3], g[3];
  v0 -= p;
  v1 -= p;
  v2 -= p;

  crossproduct(v1 - v0, v2 - v0, g);
  if (dot(n, v0) <= 0 && dot(n, v1) <= 0 && dot(n, v2) <= 0)
    return 0;
  if (dot(g, v0) >= 0)
    return 0;

  normalise(v0);
  normalise(v1);
  normalise(v2);

  crossproduct(v1 - v0, v2 - v0, a);
  return MAX(-dotproduct(n, a), 0.f);
}


static inline mf_t
lh_sample_light(
    path_t *path,
    uint64_t prim_beg,
    int prim_cnt,
    const float *p,
    const float *n)
{
  // pdf /= num_lights;
  // return offset + clamp(int(rng * num_lights), 0, int(num_lights - 1));

  const int stride = ceilf(num_lights, LH_MAX_EVALS);
  const float rng = pointsampler(path, s_dim_nee_light2);

  int lane;
  {
    float lane_ = MIN(floorf(rng * stride), stride - 1);
    lane = (int)lane_;
    rng -= lane_ / stride;
    rng *= stride;
    pdf /= stride;
  }

  int selected_idx = -1;
  float selected_mass = 0;
  float mass = 0;

  for (int i = offset + lane; i < offset + num_lights; i += stride)
  {
    lh_prim_t *lhprim = rt.lights->prims + i;
    float m = lh_projected_tri_area(lhprim, p, n);

    if (m > 0)
    {
      float f = mass;
      mass += m;
      f /= mass;

      if (rng < f)
      {
        rng /= f;
      }
      else
      {
        selected_idx = i;
        selected_mass = m;
        rng = (rng - f) / (1. - f);
      }
      rng = CLAMP(rng, 0., 1.);
    }
  }

  mass /= pdf;
  pdf = selected_mass;

  if (!(pdf > 0))
      return mf_set1(0.0f);

  pdf /= mass;


  // TODO: go to primid and sample primitive via
  float rx = pointsampler(path, s_dim_nee_x);
  float ry = pointsampler(path, s_dim_nee_y);

  // TODO: fill vertex on path
  int v = p->length;

  // sample point on triangle
  p->v[v].hit.prim = rt.lights->prim[sampled_idx].primid;
  prims_sample(rt.prims, rt.lights->primid[sampled_idx], rx, ry, &path->v[v].hit, path->time);

  if(v)
  {
    // only init segment if we're not starting a light ray
    for(int k=0;k<3;k++)
      path->e[v].omega[k] = path->v[v].hit.x[k] - path->v[v-1].hit.x[k];
    path->e[v].dist = sqrtf(dotproduct(path->e[v].omega, path->e[v].omega));
    for(int k=0;k<3;k++) path->e[v].omega[k] *= 1./path->e[v].dist;
  }
  // init emission and normals and such
  shader_prepare(path, v);
  
  path->v[v].pdf = mf_set1(pdf);
  if(path->v[v].shading.roughness > 1.0f-1e-4f)
    path->v[v].material_modes = path->v[v].mode = s_emit | s_diffuse;
  else
    path->v[v].material_modes = path->v[v].mode = s_emit | s_glossy;
  return mf_div(path->v[v].shading.em, path->v[v].pdf);
}

static inline void
lh_aabb_vertex(float aabb[6], int i, float *out)
{
  out[0] = aabb[3*((i>>0) & 1) + 0];
  out[1] = aabb[3*((i>>1) & 1) + 1];
  out[2] = aabb[3*((i>>2) & 1) + 2];
}

static inline float
_acos(float x)
{
  //return acos(clamp(x, -0.9999999, 0.9999999));

  float a = -0.939115566365855;
  float b = 0.9217841528914573;
  float c = -1.2845906244690837;
  float d = 0.295624144969963174;
  return 3.1415926535 * 0.5 + (a * x + b * x * x * x) / (1 + c * x * x + d * x * x * x * x);

  //return 1.57079f - 1.57079f * x; /* not good enough */
}

static inline float
lh_importance_measure(
    const lh_node_t *root,
    const int c,
    const float *p,
    const float *n)
{
  // centroid of the aabb
  // TODO:
  vec3 c = 0.5f * (root.aabb[0] + root.aabb[1]);
  vec3 pc = c - p;
  vec3 pc_n = normalise(pc);
  vec3 dc2 = 0.5 * (c - root.aabb[0]);
  float r2_sq = dotproduct(dc2, dc2);
  float d_sq = MAX(r2_sq, dot(pc, pc));

  float th_s = cone.th_o + cone.th_e;
  float th = _acos(dot(cone.axis, -pc_n));
  float _max = -1.0f/0.0f;
  float th_i = _acos(dot(n, pc_n));
  float m_d = 1.0/0.0;
  vec3 m_cu;
  float cos_th_u = 1.0f;

  for(int i=0;i<8;i++)
  {
    float e[3];
    lh_aabb_vertex(curr->aabb[c], i, e);
    float pe_n[3] = {e[0]-p[0],e[1]-p[1],e[2]-p[2]};
    normalise(pe_n);
    _max = MAX(_max, dotproduct(n, pe_n));
    float d = dotproduct(e, cone.axis);
    m_d = MIN(m_d, d);
    m_cu = m_d == d ? e : m_cu;
    cos_th_u = MIN(cos_th_u, dotproduct(pc_n, pe_n));
  }

  vec3 cup = p - m_cu;
  #if LH_PLANE_CHECKS
  if(_max <= 0 || (th_s <= lh_b2v(128, 0, M_PI) &&
     dotproduct(cup, cone.axis) <= 0))
    return 0.0f;
  #endif

  float th_u = _acos(cos_th_u);

  float th_h = MAX(th - cone.th_o - th_u, 0.0f);
  float th_ih = MAX(th_i - th_u, 0.0f);

  float cos_th_h = th_h < cone.th_e ? cos(th_h) : 0.00f;
  return (fabsf(cosf(th_ih)) * root->energy_combined * cos_th_h) / d_sq;
}



// returns throughput into direction of segment, no visibility is checked here.
// will construct vertex p->v[p->length] and not increment length
mf_t lights_sample_next_event(path_t *path)
{
  // TODO: traverse light hierarchy
  // int
// lh_sample_light_lh(vec3 p, vec3 n, float rng, inout float pdf)
// {
  float rng = pointsampler(path, s_dim_nee_light1);
  uint64_t nodei = 0;
  lh_node_t *curr = rt.lights->lh.nodes;

  while (1)
  {
    float acc = 0.0f;
    float cdf[4];
    int last = 0;

    // normal n needs to be normalised
    for(int i=0;i<4;i++)
    {
      if (root->child[i])
      {
        lh_node_t root = rt.light->nodes + curr->child[i];
        acc += lh_importance_measure(curr, i, p, n);
        last = i;
      }
      cdf[i] = acc;
    }

    if (acc == 0.0f)
      return mf_set1(0.0f);

    float prev_thresh = 0.0f;
    for (int c = 0; c < 4; c++)
    {
      uint64_t child = curr->child[c];
      float next_thresh = cdf[c] / acc;
      int prim_cnt = (child & ~(1ul<<63)) >> 55; // primitive count, currently one or empty leaf
      uint64_t prim_beg = child & ((1<<55)-1);   // primitive start index
      if(rng < next_thresh || c == last)
      {
        float diff = next_thresh - prev_thresh;
        rng = (rng - prev_thresh) / diff;
        pdf *= diff;

        if(child & (1ul<<63))
        {
          return lh_sample_light(path, prim_beg, prim_cnt, p, n);
        }
        else
        {
          curr = rt.lights->lh.nodes + child;
          break;
        }
      }
      prev_thresh = next_thresh;
    }
  }
}

// sample a light ray at p->v[0] and p->e[1]
mf_t lights_sample(path_t *p)
{
  assert(0 && "TODO implement light tracing");
  return mf_set1(0.0f);
}

// returns the pdf of sampling the first or last segment
// of the path (depending on v) via lights_sample()
mf_t lights_pdf(const path_t *p, int v)
{
  assert(0 && "TODO implement light tracing");
  return mf_set1(0.0f);
}

// evaluate emitted radiance of given path vertex v,
// in direction towards the sensor.
mf_t lights_eval_vertex(path_t *path, int v)
{
  // TODO: copy from list.c or put in shared header
}
