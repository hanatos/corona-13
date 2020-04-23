/*
    This file is part of corona-13.

    copyright (c) 2015 johannes hanika.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "prims.h"
#include "accel.h"

#include <math.h>
#include <float.h>
#include <assert.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define MAX_TREE_DEPTH 40
#define SAH_TESTS 7
#define SAH_LOG_STEP 3
// has to be < 32
#define NUM_TRIS_PER_LEAF 6
// if no split is found, a leaf is made if tricount < 32.
//#define STUPIDLY_SPLIT_AT_RANDOM  // ..or partition object list in middle..
//#define KICK_SILLY_TRIANGLES  // ..or just make a leaf with the first 31 triangles ;)

// #define BVH_SSE_SAH

void accel_cleanup(accel_t *b)
{
  free(b->tree);
  if(b->shadow_cache) free(b->shadow_cache);
  free(b);
}

void checktree(accel_t *b, qbvh_node_t *node, qbvh_stats_t *stats, const int depth)
{
  if(depth == 0) stats->parent_box = (b->aabb[5] - b->aabb[2])*(b->aabb[4] - b->aabb[1]) + (b->aabb[3] - b->aabb[0])*(b->aabb[4] - b->aabb[1]) + (b->aabb[5] - b->aabb[2])*(b->aabb[3] - b->aabb[0]);
  const float pbox = stats->parent_box;
  int child_tris[4] = {0,0,0,0};
  float child_box[4];
  for(int c=0;c<4;c++)
    child_box[c] = (node->aabb0[5].f[c] - node->aabb0[2].f[c])*(node->aabb0[4].f[c] - node->aabb0[1].f[c]) +
                   (node->aabb0[3].f[c] - node->aabb0[0].f[c])*(node->aabb0[4].f[c] - node->aabb0[1].f[c]) +
                   (node->aabb0[5].f[c] - node->aabb0[2].f[c])*(node->aabb0[3].f[c] - node->aabb0[0].f[c]);
  stats->num_nodes++;
  // recurse in all children, collect num tris etc.
  assert(node->axis0  >= 0 && node->axis0  < 3);
  assert(node->axis00 >= 0 && node->axis00 < 3);
  assert(node->axis01 >= 0 && node->axis01 < 3);
  assert(node - b->tree == 0 || (node->parent != -1 && node->parent < b->num_nodes));
  for(int c=0;c<4;c++)
  {
    int tris_before = stats->num_tris;
    // leaf?
    if(node->child[c] & (0x80000000ul<<32))
    {
      stats->num_leaves++;
      stats->sum_depth += depth;
      if(depth > stats->max_depth) stats->max_depth = depth;
      const int num_prims = node->child[c] & ~-(1<<5);
      const int prims = (node->child[c] ^ (0x80000000ul<<32)) >> 5;
      if(num_prims > 0)
      {
        stats->num_tris += num_prims;
        if(prims + num_prims > (int)b->prims->num_prims) fprintf(stderr, "[accel] ERR: prims index out of range!\n");
      }
      else stats->num_empties++;
    }
    else
    {
      if(child_box[c] > 0) stats->parent_box = child_box[c];
      else                 stats->parent_box = 0.0f;
      qbvh_node_t *child = b->tree + node->child[c];
      checktree(b, child, stats, depth + 1);
    }
    child_tris[c] = stats->num_tris - tris_before;
    if(isfinite(child_box[c]) && isfinite(pbox)) stats->SA += (child_tris[c]*child_box[c])/pbox;
  }
}

static inline void accel_refit(accel_t *b, qbvh_node_t *node)
{
  for(int c=0;c<4;c++)
  {
    // leaf?
    if(node->child[c] & (0x80000000ul<<32))
    {
      for(int k=0;k<3;k++) node->aabb1[k].f[c] = FLT_MAX;
      for(int k=0;k<3;k++) node->aabb1[3+k].f[c] = -FLT_MAX;
      const uint64_t num_prims = node->child[c] & ~-(1<<5);
      const uint64_t prims = (node->child[c] ^ (0x80000000ul<<32)) >> 5;
      for(uint64_t k=prims;k<prims+num_prims;k++)
      {
        // get bounds from geometric primitive and refit shutter close box aabb1[c]
        for(int d=0;d<3;d++)
        {
          float m, M;
          prims_get_bounds_shutter_close(b->prims, b->prims->primid[k], d, &m, &M);
          node->aabb1[d].f[c] = fminf(node->aabb1[d].f[c], m);
          node->aabb1[3+d].f[c] = fmaxf(node->aabb1[3+d].f[c], M);
        }
      }
    }
    else
    {
      qbvh_node_t *child = b->tree + node->child[c];
      accel_refit(b, child);
      // set aabb1[c] to fit the four boxes of child[c]
      for(int d=0;d<6;d++)
        node->aabb1[d].f[c] = child->aabb1[d].f[0];
      for(int k=1;k<4;k++)
        for(int d=0;d<3;d++)
        {
          node->aabb1[d].f[c] = fminf(node->aabb1[d].f[c], child->aabb1[d].f[k]);
          node->aabb1[3+d].f[c] = fmaxf(node->aabb1[3+d].f[c], child->aabb1[3+d].f[k]);
        }
    }
  }
}

accel_t* accel_init(prims_t *prims)
{
  accel_t *b = (accel_t*)malloc(sizeof(accel_t));
  b->shadow_cache = NULL;
  b->aabb[0] = b->aabb[1] = b->aabb[2] = FLT_MAX;
  b->aabb[3] = b->aabb[4] = b->aabb[5] = - FLT_MAX;
  b->prims = prims;

  b->num_nodes = 0;
  b->node_bufsize = b->prims->num_prims > 10 ? b->prims->num_prims: 10;
  b->tree = (qbvh_node_t*)common_alloc(128, b->node_bufsize * sizeof(qbvh_node_t));
  if(!b->tree)
    fprintf(stderr, "qbvh could not allocate enough memory for node buffer!\n");
#ifdef BVH_SSE_SAH
  b->tri_aabb = (__m128 *)common_alloc(128, sizeof(float) * 4 * 2 * b->prims->num_prims);
  if(!b->tri_aabb)
    fprintf(stderr, "qbvh could not allocate enough memory for aabb buffer!\n");
#endif

  // init shadow cache:
  b->shadow_cache_last = (1<<12) - 1;
  b->shadow_cache = (unsigned int *)malloc(rt.num_threads * sizeof(unsigned int) * (b->shadow_cache_last+1));
  memset(b->shadow_cache, 0, rt.num_threads * sizeof(unsigned int) * (b->shadow_cache_last+1));

  return b;
}

float getSplitDim(accel_t *b, const int left_in, const int right_in, const float aabb[6],
    const int d, float *split)
{
  //*split = 0.5f*(aabb[3+d] + aabb[d]);
  //return 0.0f;

  float min, max;

  const int numc = SAH_TESTS;
  int binmin[numc+1], binmax[numc+1];
  const int p = d == 2 ? 0 : d+1;
  const int q = d == 0 ? 2 : d-1;
  float best = (aabb[d] + aabb[3+d])*0.5f;
  float best_score = FLT_MAX;
  const float width = aabb[d+3] - aabb[d];
  for(int k=0;k<numc+1;k++)
    binmin[k] = binmax[k] = 0;
  const int step = (int)(log10f(SAH_LOG_STEP*(right_in - left_in) + 1.0f) + 1.0f);
  for(int k=left_in;k<right_in;k+=step)
  {
#ifdef BVH_SSE_SAH
    min = ((float *)b->tri_aabb)[8 * b->data[k] + d];
    max = ((float *)b->tri_aabb)[8 * b->data[k] + 4 + d];
#else
    prims_get_bounds_shutter_open(b->prims, b->prims->primid[k], d, &min, &max);
#endif
    int imin = (int)((numc+1)*fmaxf(fminf(0.999f, (min - aabb[d])/width), 0.0f));
    int imax = (int)((numc+1)*fmaxf(fminf(0.999f, (max - aabb[d])/width), 0.0f));
    binmin[imin]++;
    binmax[imax]++;
  }
  // this is kd-code ;)
  int left = binmin[0];
  int right = 0;
  const float com = (aabb[3+p] - aabb[p])*(aabb[3+q] - aabb[q]);
  for(int i=1;i<numc+1;i++) right += binmax[i];
  for(int k=0;k<numc;k++)
  {
    float splitc = aabb[d] + (k+1)*width/(numc + 1.0f);
    float stepl = com +
        (aabb[3+p] - aabb[p])*(splitc - aabb[d]) +
        (aabb[3+q] - aabb[q])*(splitc - aabb[d]);
    float stepr = com +
        (aabb[3+p] - aabb[p])*(aabb[3+d] - splitc) +
        (aabb[3+q] - aabb[q])*(aabb[3+d] - splitc);
    float score = stepl*left + stepr*right;
    if(score < best_score)
    {
      best_score = score;
      best = splitc;
    }
    left += binmin[k+1];
    right -= binmax[k+1];
  }
  *split = best;
  return best_score;
}

float getSplit(accel_t *b, const int left_in, const int right_in, const float aabb[6],
    int *dim, float *split)
{
  // get split plane
  *dim = 0;
  for(int k=1;k<3;k++)
    if((aabb[k+3] - aabb[k]) > (aabb[*dim+3] - aabb[*dim])) *dim = k;
  return getSplitDim(b, left_in, right_in, aabb, *dim, split);
}

void sort_tris_no_both(accel_t *b, const int axis, const float split, const int left, const int right, unsigned int *pivot, float *aabblo, float *aabbro)
{
  int fsplit = right;
#ifdef BVH_SSE_SAH
  __m128 aabbl[2], aabbr[2];
  aabbl[0] = aabbr[0] = _mm_set_ps1(FLT_MAX);
  aabbl[1] = aabbr[1] = _mm_set_ps1(- FLT_MAX);
#else
  float *aabbl = aabblo, *aabbr = aabbro;
  for(int k=0;k<3;k++) { aabbl[k+3] = aabbr[k+3] = -FLT_MAX; aabbl[k] = aabbr[k] = FLT_MAX; }
#endif
  float pmin, pmax;
#ifndef BVH_SSE_SAH
  const int p = axis == 2 ? 0 : axis + 1;
  const int q = axis == 0 ? 2 : axis - 1;
#endif
  for(int i=left;i<fsplit;)
  {
#ifndef BVH_SSE_SAH
    prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], axis, &pmin, &pmax);
#else
    pmin = ((float *)b->tri_aabb)[8 * b->data[i] + axis];
    pmax = ((float *)b->tri_aabb)[8 * b->data[i] + 4 + axis];
#endif
    if(0.5f*(pmin + pmax) >= split)
    {
#ifdef BVH_SSE_SAH
      aabbr[0] = _mm_min_ps(aabbr[0], b->tri_aabb[2 * b->data[i]]);
      aabbr[1] = _mm_max_ps(aabbr[1], b->tri_aabb[2 * b->data[i]+1]);
#else
      if(aabbr[axis] > pmin) aabbr[axis] = pmin;
      if(aabbr[axis+3] < pmax) aabbr[axis+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], p, &pmin, &pmax);
      if(aabbr[p] > pmin) aabbr[p] = pmin;
      if(aabbr[p+3] < pmax) aabbr[p+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], q, &pmin, &pmax);
      if(aabbr[q] > pmin) aabbr[q] = pmin;
      if(aabbr[q+3] < pmax) aabbr[q+3] = pmax;
#endif
      //right, rotate bufs:
      primid_t tmp = b->prims->primid[i];
      b->prims->primid[i] = b->prims->primid[--fsplit];
      b->prims->primid[fsplit] = tmp;
    }
    else
    {
#ifdef BVH_SSE_SAH
      aabbl[0] = _mm_min_ps(aabbl[0], b->tri_aabb[2 * b->data[i]]);
      aabbl[1] = _mm_max_ps(aabbl[1], b->tri_aabb[2 * b->data[i]+1]);
#else
      //left
      if(aabbl[axis] > pmin) aabbl[axis] = pmin;
      if(aabbl[axis+3] < pmax) aabbl[axis+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], p, &pmin, &pmax);
      if(aabbl[p] > pmin) aabbl[p] = pmin;
      if(aabbl[p+3] < pmax) aabbl[p+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], q, &pmin, &pmax);
      if(aabbl[q] > pmin) aabbl[q] = pmin;
      if(aabbl[q+3] < pmax) aabbl[q+3] = pmax;
#endif
      i++;
    }
  }
#ifdef BVH_SSE_SAH
  for(int i=0;i<3;i++) aabblo[i] = ((float*)aabbl)[i];
  for(int i=0;i<3;i++) aabblo[3+i] = ((float*)aabbl)[i+4];
  for(int i=0;i<3;i++) aabbro[i] = ((float*)aabbr)[i];
  for(int i=0;i<3;i++) aabbro[3+i] = ((float*)aabbr)[i+4];
#endif
  *pivot = fsplit;
}

void sort_tris(accel_t *b, const int axis, const float split, const int left, const int right, unsigned int *pivot, float *aabblo, float *aabbro)
{
  //sort_tris_no_both(b, axis, split, left, right, pivot, aabblo, aabbro);
  //return;

  int both = right;
  int fsplit = right;
#ifdef BVH_SSE_SAH
  __m128 aabbl[2], aabbr[2], aabbb[2];
  aabbb[0] = aabbl[0] = aabbr[0] = _mm_set_ps1(FLT_MAX);
  aabbb[1] = aabbl[1] = aabbr[1] = _mm_set_ps1(- FLT_MAX);
#else
  float *aabbl = aabblo, *aabbr = aabbro;
  float aabbb[6];
  for(int k=0;k<3;k++) { aabbb[k+3] = aabbl[k+3] = aabbr[k+3] = -FLT_MAX; aabbb[k] = aabbl[k] = aabbr[k] = FLT_MAX; }
#endif
  float pmin, pmax;
#ifndef BVH_SSE_SAH
  const int p = axis == 2 ? 0 : axis + 1;
  const int q = axis == 0 ? 2 : axis - 1;
#endif
  for(int i=left;i<both;)
  {
#ifndef BVH_SSE_SAH
    prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], axis, &pmin, &pmax);
#else
    pmin = ((float *)b->tri_aabb)[8 * b->data[i] + axis];
    pmax = ((float *)b->tri_aabb)[8 * b->data[i] + 4 + axis];
#endif
    if(pmin >= split)
    {
#ifdef BVH_SSE_SAH
      aabbr[0] = _mm_min_ps(aabbr[0], b->tri_aabb[2 * b->data[i]]);
      aabbr[1] = _mm_max_ps(aabbr[1], b->tri_aabb[2 * b->data[i]+1]);
#else
      if(aabbr[axis] > pmin) aabbr[axis] = pmin;
      if(aabbr[axis+3] < pmax) aabbr[axis+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], p, &pmin, &pmax);
      if(aabbr[p] > pmin) aabbr[p] = pmin;
      if(aabbr[p+3] < pmax) aabbr[p+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], q, &pmin, &pmax);
      if(aabbr[q] > pmin) aabbr[q] = pmin;
      if(aabbr[q+3] < pmax) aabbr[q+3] = pmax;
#endif
      //right, rotate bufs:
      primid_t tmp = b->prims->primid[i];
      b->prims->primid[i] = b->prims->primid[--both];
      b->prims->primid[both] = b->prims->primid[--fsplit];
      b->prims->primid[fsplit] = tmp;
    }
    else if(pmax <= split)
    {
#ifdef BVH_SSE_SAH
      aabbl[0] = _mm_min_ps(aabbl[0], b->tri_aabb[2 * b->data[i]]);
      aabbl[1] = _mm_max_ps(aabbl[1], b->tri_aabb[2 * b->data[i]+1]);
#else
      //left
      if(aabbl[axis] > pmin) aabbl[axis] = pmin;
      if(aabbl[axis+3] < pmax) aabbl[axis+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], p, &pmin, &pmax);
      if(aabbl[p] > pmin) aabbl[p] = pmin;
      if(aabbl[p+3] < pmax) aabbl[p+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], q, &pmin, &pmax);
      if(aabbl[q] > pmin) aabbl[q] = pmin;
      if(aabbl[q+3] < pmax) aabbl[q+3] = pmax;
#endif
      i++;
    }
    else
    {
#ifdef BVH_SSE_SAH
      aabbb[0] = _mm_min_ps(aabbb[0], b->tri_aabb[2 * b->data[i]]);
      aabbb[1] = _mm_max_ps(aabbb[1], b->tri_aabb[2 * b->data[i]+1]);
#else
      if(aabbb[axis] > pmin) aabbb[axis] = pmin;
      if(aabbb[axis+3] < pmax) aabbb[axis+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], p, &pmin, &pmax);
      if(aabbb[p] > pmin) aabbb[p] = pmin;
      if(aabbb[p+3] < pmax) aabbb[p+3] = pmax;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], q, &pmin, &pmax);
      if(aabbb[q] > pmin) aabbb[q] = pmin;
      if(aabbb[q+3] < pmax) aabbb[q+3] = pmax;
#endif
      //both
      primid_t tmp = b->prims->primid[i];
      b->prims->primid[i] = b->prims->primid[--both];
      b->prims->primid[both] = tmp;
    }
  }
#ifdef BVH_SSE_SAH
  if((fsplit == right) || ((both != left) && (0.5f*(((float *)aabbb)[axis] + ((float*)aabbb)[axis+4]) > split)))
  {
    *pivot = both;
    aabbr[0] = _mm_min_ps(aabbr[0], aabbb[0]);
    aabbr[1] = _mm_max_ps(aabbr[1], aabbb[1]);
  }
  else
  {
    *pivot = fsplit;
    aabbl[0] = _mm_min_ps(aabbl[0], aabbb[0]);
    aabbl[1] = _mm_max_ps(aabbl[1], aabbb[1]);
  }
  for(int i=0;i<3;i++) aabblo[i]   = ((float*)aabbl)[i];
  for(int i=0;i<3;i++) aabblo[3+i] = ((float*)aabbl)[i+4];
  for(int i=0;i<3;i++) aabbro[i]   = ((float*)aabbr)[i];
  for(int i=0;i<3;i++) aabbro[3+i] = ((float*)aabbr)[i+4];
#else
  //printf("sort: [%d %d %d ]%d\n", left, both, fsplit, right);
  // decide on B->L or B->R
  if((fsplit == right) || ((both != left) && (0.5f*(aabbb[axis] + aabbb[axis+3]) > split)))
  {
    *pivot = both;
    for(int k=0;k<3;k++)
    {
      if(aabbro[k] > aabbb[k]) aabbro[k] = aabbb[k];
      if(aabbro[k+3] < aabbb[k+3]) aabbro[k+3] = aabbb[k+3];
    }
  }
  else
  {
    *pivot = fsplit;
    for(int k=0;k<3;k++)
    {
      if(aabblo[k] > aabbb[k]) aabblo[k] = aabbb[k];
      if(aabblo[k+3] < aabbb[k+3]) aabblo[k+3] = aabbb[k+3];
    }
  }
#endif
}

void build_help(accel_t *b, qbvh_node_t *node, const int left, const int right, float *paabb, const int depth, qbvh_node_t *parent, const int child)
{
  node->parent = parent - b->tree;
  int retry = 0;
  //printf("building node %d at depth %d with [%d..]%d\n", node - b->tree, depth, left, right);
  // split axis0
  int axis0, axis00, axis01;
  unsigned int part[5];
  float split0, split1l, split1r;
  float aabb[4][6];
  part[0] = left;
  part[4] = right;
#if 0
  getSplit(b, left, right, aabb, &axis0, &split0);
#else
  axis0 = 0;
  float score, split;
  float best_score = getSplitDim(b, left, right, paabb, 0, &split0);
  if((score = getSplitDim(b, left, right, paabb, 1, &split)) < best_score) { best_score = score; split0 = split; axis0 = 1; }
  if((score = getSplitDim(b, left, right, paabb, 2, &split)) < best_score) { best_score = score; split0 = split; axis0 = 2; }

#if 1
  if(best_score < 0.0f && (right - left) < 32)
  {
    // terminate
    // FIXME: wastes mem for already alloc'ed node :(
    if(parent)
    {
      //printf("making parent a leaf.\n");
      parent->child[child] = (0x80000000ul<<32) | (left<<5) | ((right - left) & ~-(1<<5));
      return;
    }
  }
#endif

do_retry:
#endif
  //if((right - left) && !(split0 >= paabb[axis0]) && (split0 <= paabb[axis0+3])) printf("silly split plane cand: [%f | %f %f]\n", paabb[axis0], split0, paabb[3+axis0]);

  // split triangle list at axis0/split0.
  if(retry) sort_tris_no_both(b, axis0, split0, left, right, &part[2], aabb[0], aabb[1]);
  else sort_tris(b, axis0, split0, left, right, &part[2], aabb[0], aabb[1]);

  // split axis1 (2x)
#if 0
  getSplit(b, left, pivot, aabbl, &axis00, &split1l);
  getSplit(b, pivot, right, aabbr, &axis01, &split1r);
#else
  axis00 = 0;
  best_score = getSplitDim(b, left, part[2], aabb[0], 0, &split1l);
  if((score = getSplitDim(b, left, part[2], aabb[0], 1, &split)) < best_score) { best_score = score; split1l = split; axis00 = 1; }
  if((score = getSplitDim(b, left, part[2], aabb[0], 2, &split)) < best_score) { split1l = split; axis00 = 2; }
  axis01 = 0;
  best_score = getSplitDim(b, part[2], right, aabb[1], 0, &split1r);
  if((score = getSplitDim(b, part[2], right, aabb[1], 1, &split)) < best_score) { best_score = score; split1r = split; axis01 = 1; }
  if((score = getSplitDim(b, part[2], right, aabb[1], 2, &split)) < best_score) { split1r = split; axis01 = 2; }
#endif

  //if((part[2] - left) && !(split1l >= paabb[axis00]) && (split1l <= paabb[axis00+3])) printf("silly split plane cand: [%f | %g %f]\n", paabb[axis00], split1l, paabb[3+axis00]);
  if(retry) sort_tris_no_both(b, axis00, split1l, left, part[2], &part[1], aabb[0], aabb[1]);
  else sort_tris(b, axis00, split1l, left, part[2], &part[1], aabb[0], aabb[1]);

  //if((right - part[2]) && !(split1r >= paabb[axis01]) && (split1r <= paabb[axis01+3])) printf("silly split plane cand: [%f | %g %f]\n", paabb[axis01], split1r, paabb[3+axis01]);
  if(retry) sort_tris_no_both(b, axis01, split1r, part[2], right, &part[3], aabb[2], aabb[3]);
  else sort_tris(b, axis01, split1r, part[2], right, &part[3], aabb[2], aabb[3]);

  if((right   - part[3] == (unsigned int)(right - left)) ||
     (part[3] - part[2] == (unsigned int)(right - left)) ||
     (part[2] - part[1] == (unsigned int)(right - left)) ||
     (part[1] - left    == (unsigned int)(right - left)))
  {
    if(!retry)
    {
      //printf("restart!\n");
      retry = 1;
      goto do_retry;
    }
#ifdef STUPIDLY_SPLIT_AT_RANDOM
    if(1)
#else
#ifdef KICK_SILLY_TRIANGLES
    if(0)
#else
    if((parent && (right - left > ~-(1<<5))) || !parent)
#endif
#endif
    {
      // printf("splitting at random: [%d %d %d %d] %d\n", part[0], part[1], part[2], part[3], part[4]);
      // split list in the middle, for 1tri/leaf case:
      part[2] = (right   + left)/2;
      part[3] = (right   + part[2])/2;
      part[1] = (part[2] + left)/2;
      // find aabbs
      for(int k=0;k<3;k++) for(int i=0;i<4;i++) { aabb[i][k] = FLT_MAX; aabb[i][k+3] = -FLT_MAX; }
      for(int p=0;p<4;p++)
      {
        for(unsigned int k=part[p];k<part[p+1];k++)
        {
          float min, max;
          for(int i=0;i<3;i++)
          {
            prims_get_bounds_shutter_open(b->prims, b->prims->primid[k], i, &min, &max);
            if(min < aabb[p][i  ]) aabb[p][i  ] = min;
            if(max > aabb[p][i+3]) aabb[p][i+3] = max;
          }
        }
      }
    }
    else
    {
      //printf("making parent a leaf.\n");
      parent->child[child] = (0x80000000lu<<32) | (left<<5) | ((right - left) & ~-(1<<5));
      return;
    }
  }

  node->axis0  = axis0;
  node->axis00 = axis00;
  node->axis01 = axis01;
  for(int k=0;k<6;k++) for(int p=0;p<4;p++) node->aabb0[k].f[p] = aabb[p][k];

  // this will at times segfault because it runs over the node array buffer.
// #define BUILD_IN_PARALLEL
// TODO: benchmark when this actually makes sense to run in parallel
#ifdef BUILD_IN_PARALLEL
  const int num_nodes = b->num_nodes;
  b->num_nodes += 4;
#pragma omp parallel for if(right-left > 20000)
#endif
  for (int p=0;p<4;p++)
  {
    // TODO: quads don't count as one but as two primitives here!
    if (depth == MAX_TREE_DEPTH || part[p+1] - part[p] <= NUM_TRIS_PER_LEAF || b->num_nodes >= b->node_bufsize - 1)
      node->child[p] = (0x80000000lu<<32) | (part[p]<<5) | ((part[p+1] - part[p]) & ~-(1<<5));
    else
    {
#ifdef BUILD_IN_PARALLEL
      node->child[p] = num_nodes + p;
#else
      node->child[p] = b->num_nodes ++;
#endif
      assert(b->num_nodes < b->node_bufsize);
      build_help(b, b->tree + node->child[p], part[p], part[p+1], aabb[p], depth + 1, node, p);
    }
  }
}

void accel_build(accel_t *b, const char* filename)
{
  b->num_nodes = 1;
  b->aabb[0] = b->aabb[1] = b->aabb[2] = FLT_MAX;
  b->aabb[3] = b->aabb[4] = b->aabb[5] = - FLT_MAX;
  for(unsigned int i=0;i<b->prims->num_prims;i++)
  {
    for(int k=0;k<3;k++)
    {
      float min, max;
      prims_get_bounds_shutter_open(b->prims, b->prims->primid[i], k, &min, &max);
#ifdef BVH_SSE_SAH
      ((float *)b->tri_aabb)[8*i + k] = min;
      ((float *)b->tri_aabb)[8*i + k + 4] = max;
#endif
      if(b->aabb[k] > min) b->aabb[k] = min;
      if(b->aabb[k+3] < max) b->aabb[k+3] = max;
    }
  }
  build_help(b, b->tree, 0, b->prims->num_prims, b->aabb, 0, NULL, 0);
  b->tree[0].parent = -1;
#ifdef BVH_SSE_SAH
  free(b->tri_aabb);
#endif
  accel_refit(b, b->tree);
#if 1
  qbvh_stats_t stats;
  memset(&stats, 0, sizeof(qbvh_stats_t));
  checktree(b, b->tree, &stats, 0);
  printf("[accel] num nodes created: %d\n", b->num_nodes);
  printf("[accel] created qbvh with %d nodes, %d leaves (avg %.3f/%.3f depth/num prims, %d empty), %d prims, max depth is %d.\n",
      stats.num_nodes, stats.num_leaves, stats.sum_depth/(float)stats.num_leaves,
      stats.num_tris/(float)(stats.num_leaves-stats.num_empties), stats.num_empties, stats.num_tris, stats.max_depth);
  printf("[accel] overall SAV: %f\n", stats.SA);
#endif
}

static inline void aabb_intersect(
    const qbvh_float4_t *aabb0,
    const qbvh_float4_t *aabb1,
    const __m128 *t0,
    const __m128 *t1,
    const __m128 *pos4,
    const __m128 *invdir4,
    const int *near,
    const int *far,
    qbvh_float4_t *tmin,
    qbvh_float4_t *tmax)
{
  for(int k=0;k<3;k++)
  {
    __m128 t[2];
    t[0] = _mm_mul_ps(_mm_sub_ps(_mm_add_ps(
            _mm_mul_ps(aabb0[k+3*near[k]].m, *t0),
            _mm_mul_ps(aabb1[k+3*near[k]].m, *t1)), pos4[k]), invdir4[k]);
    t[1] = _mm_mul_ps(_mm_sub_ps(_mm_add_ps(
            _mm_mul_ps(aabb0[k+3*far [k]].m, *t0),
            _mm_mul_ps(aabb1[k+3*far [k]].m, *t1)), pos4[k]), invdir4[k]);
    tmin->m = _mm_max_ps(tmin->m, t[0]);
    tmax->m = _mm_min_ps(tmax->m, t[1]);
  }
}

void accel_intersect(const accel_t *b, const ray_t *ray, hit_t *hit)
{
  // ((accel_t *)b)->boxtests++;
  int near[3];
  int far[3];
  for(int k=0;k<3;k++)
  {
    // need to take sign bit or test invdir for < 0 to catch +-inf cases in aabb intersection.
    near[k] = (*(unsigned int*)(&(ray->dir[k])))>>31;
    far[k] = 1 ^ near[k];
  }

#if 0
  int dir = 0;
  if(fabsf(ray->dir[1]) > fabsf(ray->dir[0])) dir = 1;
  if(fabsf(ray->dir[2]) > fabsf(ray->dir[dir])) dir = 2;
  if(ray->dir[dir] < 0.0f) dir += 3;
  dir *= 5;
#endif

  const qbvh_node_t *node = b->tree;
  uint64_t stack[3*MAX_TREE_DEPTH];
  int stackpos = 0;
  uint64_t current;
  stack[0] = 0;

  __m128 invdir4[3], pos4[3];
  const __m128 t0 = _mm_set1_ps(1.0f-ray->time);
  const __m128 t1 = _mm_set1_ps(ray->time);

  for(int k=0;k<3;k++)
  {
    invdir4[k] = _mm_set1_ps(1.0f/ray->dir[k]);
    pos4[k] = _mm_load_ps1(ray->pos + k);
  }

  while(1)
  {
    // 4x aabb test, sort and push selected index to stack, hold one reg
    // push children according to volume intersection test.
    qbvh_float4_t tmin4, tmax4;
    tmin4.m = _mm_setzero_ps();
    tmax4.m = _mm_set_ps1(hit->dist);
    //((accel_t *)b)->boxtests++;
    aabb_intersect(node->aabb0, node->aabb1, &t0, &t1, pos4, invdir4, near, far, &tmin4, &tmax4);
    const __m128 leq = _mm_cmple_ps(tmin4.m, tmax4.m);
    const unsigned int *i = (unsigned int*)&leq;
    //const unsigned int i[4] = {1, 1, 1, 1};
    const int axis0  = node->axis0;
    const int axis1n = near[axis0] ? node->axis01 : node->axis00;
    const int axis1f = near[axis0] ? node->axis00 : node->axis01;

    const unsigned int n11 = (far [axis0]<<1) | far [axis1f];
    const unsigned int n10 = (far [axis0]<<1) | near[axis1f];
    const unsigned int n01 = (near[axis0]<<1) | far [axis1n];
    const unsigned int n00 = (near[axis0]<<1) | near[axis1n];
   
    if(i[n00])
    {
      current = node->child[n00];
      if(i[n11]) stack[stackpos++] = node->child[n11];
      if(i[n10]) stack[stackpos++] = node->child[n10];
      if(i[n01]) stack[stackpos++] = node->child[n01];
    }
    else if(i[n01])
    {
      current = node->child[n01];
      if(i[n11]) stack[stackpos++] = node->child[n11];
      if(i[n10]) stack[stackpos++] = node->child[n10];
    }
    else if(i[n10])
    {
      current = node->child[n10];
      if(i[n11]) stack[stackpos++] = node->child[n11];
    }
    else if(i[n11]) current = node->child[n11];
    else
    {
      if(stackpos == 0) return;
      stackpos--;
      current = stack[stackpos];
    }

    while(current & (0x80000000ul<<32))
    {
      // intersect primitives
      uint64_t tmp_index = (current ^ (0x80000000ul<<32)) >> 5;
      const uint64_t num = current & ~-(1<<5);
      for (uint64_t i = 0; i < num; i++)
      {
        prims_intersect(b->prims, b->prims->primid[tmp_index], ray, hit);
        tmp_index++;
      }
      if (stackpos == 0) return;
      --stackpos;
      current = stack[stackpos];
    }
    node = b->tree + current;
  }
}

int accel_visible(const accel_t *b, const ray_t *ray)
{
  //((accel_t *)b)->boxtests++;
  int near[3];
  int far[3];
  for(int k=0;k<3;k++)
  {
    near[k] = (*(unsigned int*)(&(ray->dir[k])))>>31;
    far[k] = 1 ^ near[k];
  }

  const int hash = (b->shadow_cache_last+1)*common_get_threadid() + (((((int)ray->pos[0])*73856093) ^ (((int)ray->pos[1])*19349663) ^ (((int)ray->pos[2])*83492791)) & b->shadow_cache_last);
  const qbvh_node_t *node = b->tree + b->shadow_cache[hash];
  const qbvh_node_t *ep = node;
  uint64_t stack[3*MAX_TREE_DEPTH];
  int stackpos = 0;
  uint64_t current;
  stack[0] = 0;

  __m128 invdir4[3], pos4[3];
  const __m128 t0 = _mm_set1_ps(1.0f-ray->time);
  const __m128 t1 = _mm_set1_ps(ray->time);

  for(int k=0;k<3;k++)
  {
    invdir4[k] = _mm_set1_ps(1.0f/ray->dir[k]);
    pos4[k] = _mm_load_ps1(ray->pos + k);
  }

  while(1)
  {
    qbvh_float4_t tmin4, tmax4;
    tmin4.m = _mm_setzero_ps();
    tmax4.m = _mm_set_ps1(1.0f);
    aabb_intersect(node->aabb0, node->aabb1, &t0, &t1, pos4, invdir4, near, far, &tmin4, &tmax4);
    const __m128 leq = _mm_cmple_ps(tmin4.m, tmax4.m);
    const int *i = (int*)&leq;
    for(int k=0;k<4;k++)
    {
      stack[stackpos] = node->child[k];
      stackpos -= i[k];
    }
    while(1)
    {
      while (stackpos == 0)
      {
        // push children of node->parent
        if(ep->parent == (uint64_t)-1) return 1;
        node = b->tree + ep->parent;
        // only push intersecting aabbs.
        qbvh_float4_t tmin4, tmax4;
        tmin4.m = _mm_setzero_ps();
        tmax4.m = _mm_set_ps1(1.0f);
        aabb_intersect(node->aabb0, node->aabb1, &t0, &t1, pos4, invdir4, near, far, &tmin4, &tmax4);
        const __m128 leq = _mm_cmple_ps(tmin4.m, tmax4.m);
        const int *i = (int*)&leq;
        for(int k=0;k<4;k++)
        {
          if(node->child[k] != ep - b->tree)
          {
            stack[stackpos] = node->child[k];
            stackpos -= i[k];
          }
        }
        ep = node;
      }
      --stackpos;
      current = stack[stackpos];
      if(!(current & (0x80000000ul<<32))) break;
      // intersect primitives
      uint64_t tmp_index = (current ^ (0x80000000ul<<32)) >> 5;
      const uint64_t num = current & ~-(1<<5);
      for (uint64_t i = 0; i < num; i++)
      {
        if(prims_intersect_visible(b->prims, b->prims->primid[tmp_index], ray))
        {
          // cache node:
          b->shadow_cache[hash] = node - b->tree;
          return 0;
        }
        tmp_index++;
      }
    }
    node = b->tree + current;
  }
}

