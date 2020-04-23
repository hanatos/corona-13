/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "corona_common.h"
#include "prims.h"
#include "accel.h"

#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


static inline float min(const float a, const float b) { return a < b ? a : b; }
static inline float max(const float a, const float b) { return a < b ? b : a; }

accel_t* accel_init(prims_t *p)
{
  accel_t *b = (accel_t *)malloc(sizeof(accel_t));
  b->aabb[0] = b->aabb[1] = b->aabb[2] = INFINITY;
  b->aabb[3] = b->aabb[4] = b->aabb[5] = - INFINITY;
  b->num_nodes = 0;
  b->prims = p;

  b->node_bufsize = b->prims->num_prims + 2;
  b->tree = (bih_node_t *)malloc(b->node_bufsize * sizeof(bih_node_t));

  return b;
}


void accel_cleanup(accel_t *b)
{
  free(b->tree);
}

int build_help(accel_t *b, bih_node_t *root, const unsigned int left, const unsigned int right, float aabb[6], float sbb[6], const int depth)
{

  float split, clipl, clipr, max, pmin, pmax, amin, amax;
  int dim = 0, tree_depth = depth;
  unsigned int fsplit;

  const int max_tree_depth = 100;


  int retry = 0;
try_harder:
  if(retry > 7)
  {
    root->data = (left<<3) | 3;
    *(unsigned int*)root = right - left;
    return depth;
  }

  for(int i=0;i<3;i++)
  {
    max = sbb[3] - sbb[0];
    dim = 0;
    for(int k=1;k<3;k++)
    {
      if(sbb[3+k] - sbb[k] > max)
      {
        max = sbb[3+k] - sbb[k];
        dim = k;
      }
    }
    split = .5f*(sbb[3+dim] + sbb[dim]);
    if(split < aabb[dim]) sbb[dim] = split;
    else if(split > aabb[3+dim]) sbb[dim+3] = split;
    else break;
  }
  fsplit = right;
  unsigned int both = right;
  float bothM = clipl = amax = -INFINITY;
  float bothm = clipr = amin = INFINITY;
  for(unsigned int i=left;i<both;)
  {
    prims_get_bounds(b->prims, b->prims->primid[i], dim, &pmin, &pmax);
    if(pmin < amin) amin = pmin;
    if(pmax > amax) amax = pmax;
    if(pmin >= split)
    {
      primid_t tmp = b->prims->primid[i];
      b->prims->primid[i] = b->prims->primid[--both];
      b->prims->primid[both] = b->prims->primid[--fsplit];
      b->prims->primid[fsplit] = tmp;
      if(pmin < clipr) clipr = pmin;
    }
    else if(pmax <= split)
    {
      if(pmax > clipl) clipl = pmax;
      i++;
    }
    else
    {
      primid_t tmp = b->prims->primid[i];
      b->prims->primid[i] = b->prims->primid[--both];
      b->prims->primid[both] = tmp;
      if(pmin < bothm) bothm = pmin;
      if(pmax > bothM) bothM = pmax;
    }
  }

  if(0.5f*(bothm + bothM) > split)
  {
    fsplit = both;
    if(bothm < clipr) clipr = bothm;
  }
  else
  {
    if(bothM > clipl) clipl = bothM;
  }
  if(1.3f*(amax - amin) < aabb[dim+3] - aabb[dim])
  {
    root->data = (b->num_nodes << 3) | 4 | dim;
    root->clip[0] = amin;
    root->clip[1] = amax;
    assert(b->node_bufsize > b->num_nodes);
    b->num_nodes++;
    b->num_clipnodes++;
    float aabb_new[6];
    for(int k=0;k<6;k++) aabb_new[k] = aabb[k];
    aabb_new[dim] = amin;
    aabb_new[dim+3] = amax;
    return build_help(b, b->tree + b->num_nodes-1, left, right, aabb_new, sbb, depth + 1);
  }
  if(fsplit == left && clipl == -INFINITY)
  {
    sbb[dim] = split;
    retry++;
    goto try_harder;
  }
  if(fsplit == right && clipr == INFINITY)
  {
    sbb[dim+3] = split;
    retry++;
    goto try_harder;
  }

  root->clip[0] = clipl;
  root->clip[1] = clipr;
  root->data = (b->num_nodes << 3) | dim;
  assert(b->node_bufsize >= b->num_nodes + 2);
  int child = b->num_nodes;
  b->num_nodes += 2;

  float aabb_new[6], sbb_new[6];
  const unsigned int tri_per_leaf = 6;

  if(right - fsplit == 0)
  {
    root->clip[1] = INFINITY;
    b->num_nodes--;
  }

  if(fsplit - left == 0)
  {
    root->clip[0] = - INFINITY;
    root->data -= 8;
    child--;
    b->num_nodes--;
  }
  else if(fsplit - left <= tri_per_leaf || depth >= max_tree_depth - 1)
  {
    b->tree[child].data = (left<<3) | 3;
    *(unsigned int*)(b->tree + child) = fsplit - left;
  }
  else
  {
    for(int k=0;k<6;k++) { aabb_new[k] = aabb[k]; sbb_new[k] = sbb[k]; }
    aabb_new[3+dim] = clipl;
    sbb_new[3+dim] = split;
    tree_depth = build_help(b, b->tree + child, left, fsplit, aabb_new, sbb_new, depth + 1);
  }

  if(right - fsplit > 0)
  {
    if(right - fsplit <= tri_per_leaf || depth >= max_tree_depth - 1)
    {
      b->tree[child + 1].data = (fsplit<<3) | 3;
      *(unsigned int*)(b->tree + child + 1) = right - fsplit;
    }
    else
    {
      for(int k=0;k<6;k++) { aabb_new[k] = aabb[k]; sbb_new[k] = sbb[k]; }
      aabb_new[dim] = clipr;
      sbb_new[dim] = split;
      int d = build_help(b, b->tree + child + 1, fsplit, right, aabb_new, sbb_new, depth + 1);
      tree_depth = tree_depth > d ? tree_depth : d;
    }
  }
  return tree_depth;
}


void accel_build(accel_t *b, const char *filename)
{
  b->num_nodes = 1;
  b->num_clipnodes = 0;

  b->aabb[0] = b->aabb[1] = b->aabb[2] = INFINITY;
  b->aabb[3] = b->aabb[4] = b->aabb[5] = - INFINITY;
  for(unsigned int i=0;i<b->prims->num_prims;i++)
  {
    for(int k=0;k<3;k++)
    {
      float min, max;
      prims_get_bounds(b->prims, b->prims->primid[i], k, &min, &max);
      if(b->aabb[k] > min) b->aabb[k] = min;
      if(b->aabb[k+3] < max) b->aabb[k+3] = max;
    }
  }
  float sbb[6];
  memcpy(sbb, b->aabb, 6*sizeof(float));

  build_help(b, b->tree, 0U, b->prims->num_prims, b->aabb, sbb, 0);
  //for(unsigned int k=0;k<b->prims->num_prims;k++) prim_init_accel(b->prims, k);
}

char aabb_intersect(const float * const aabb, const ray_t *ray, float *invdir, float *tmin, float *tmax)
{
  for(int k=0;k<3;k++)
  {
    float t1 = (aabb[k] - ray->pos[k]) * invdir[k];
    float t2 = (aabb[k+3] - ray->pos[k]) * invdir[k];
    if(t1 > t2)
    {
      *tmin = max(*tmin, t2);
      *tmax = min(*tmax, t1);
    }
    else
    {
      *tmin = max(*tmin, t1);
      *tmax = min(*tmax, t2);
    }
  }
  return *tmin <= *tmax;
}

#if 0
void checktree(accel_t *tree, bih_node_t *n, const int i)
{
  if(accel_is_leaf(n + 1))
    printf("number of primitives in leafnode %d: %d\n", i, *(unsigned int*)n+i);
  else
  {
    if(n[i].clip[0] != - INFINITY)

      checktree(tree, n, ((int)((((char*)bih_children(n + i, tree))) - (char*)n))/sizeof(bih_node));
    if(n[i].clip[1] != INFINITY)
      checktree(tree, n, ((int)((((char*)bih_children(n + i, tree)) + sizeof(bih_node)) - (char*)n))/sizeof(bih_node));
  }
}
#endif

void accel_intersect(const accel_t *b, const ray_t *ray, hit_t *hit)
{
  ((accel_t *)b)->boxtests++;
  const int max_tree_depth = 100;

  int below[3];
  int far[3];
  float invdir[3];
  for(int k=0;k<3;k++)
  {
    below[k] = (*(unsigned int*)(&(ray->dir[k])))>>31;
    far[k] = 1 ^ below[k];
    invdir[k] = 1.0f/ray->dir[k];
  }

  float tmin = 0.0f;
  float tmax = hit->dist;
  bih_node_t *node = b->tree;
  bih_node_t *fbSide[2];

  bih_node_t *nList[max_tree_depth];
  int nListPos = 0;
  float tminlist[max_tree_depth];
  float tmaxlist[max_tree_depth];
  int axis;
  float dist0, dist1;

  if (!aabb_intersect(b->aabb, ray, invdir, &tmin, &tmax)) return;

  while (1)
  {
    while (1)
    {
      axis = node->data & 3;
      if(axis == 3) break;
      dist0 = (node->clip[below[axis]] - ray->pos[axis]) * invdir[axis];
      dist1 = (node->clip[far[axis]] - ray->pos[axis]) * invdir[axis];
      fbSide[0] = fbSide[1] = b->tree + (node->data >> 3);//bih_children(node, b);
      fbSide[far[axis]]++;

      if(node->data & 4)//bih_isClipNode(node))
      {
        tmin = max(tmin, dist0);
        tmax = min(tmax, dist1);
        if (tmin > tmax)
        {
          if ( nListPos > 0 )
          {
            --nListPos;
            tmin = tminlist[nListPos];
            tmax = min(tmaxlist[nListPos], hit->dist);
            node = nList[nListPos];
          }
          else return;
        }
        else node = b->tree + (node->data >> 3);//bih_children(node, b);
        continue;
      }

      if ( dist0 < tmin )
      {
        if ( dist1 > tmax )
        {
          if (nListPos > 0)
          {
            --nListPos;
            tmin = tminlist[nListPos];
            tmax = min(tmaxlist[nListPos], hit->dist);
            node = nList[nListPos];
          }
          else return;
        }
        else
        {
          tmin = max(tmin, dist1);
          node = fbSide[1];
        }
      }
      else if ( dist1 > tmax )
      {
        tmax = min(tmax, dist0);
        node = fbSide[0];
      }
      else
      {
        tminlist[nListPos] = max(tmin, dist1);
        tmaxlist[nListPos] = tmax;
        nList[nListPos] = fbSide[1];
        ++nListPos;

        node = fbSide[0];
        tmax = min(tmax, dist0);
      }
    }

    unsigned int tmp_index = node->data >> 3;//bih_prims(node);
    for (unsigned int i = 0; i < *(unsigned int*)node; i++)//bih_num(node); i++)
    {
      prims_intersect(b->prims, b->prims->primid[tmp_index], ray, hit);
      tmp_index++;
    }

    while(1)
    {
      if (nListPos == 0) return;
      --nListPos;
      tmin = tminlist[nListPos];
      if (hit->dist >= tmin)
      {
        tmax = min(tmaxlist[nListPos], hit->dist);
        node = nList[nListPos];
        break;
      }
    }
  }
}

int  accel_visible(const struct accel_t *b, const ray_t *ray)
{ 
  ((accel_t *)b)->boxtests++;
  const int max_tree_depth = 100;

  int below[3];
  int far[3];
  float invdir[3];
  for(int k=0;k<3;k++)
  {
    below[k] = (*(unsigned int*)(&(ray->dir[k])))>>31;
    far[k] = 1 ^ below[k];
    invdir[k] = 1.0f/ray->dir[k];
  }

  float tmin = 0.0f;
  float tmax = 1.0f;
  bih_node_t *node = b->tree;
  bih_node_t *fbSide[2];

  bih_node_t *nList[max_tree_depth];
  int nListPos = 0;
  float tminlist[max_tree_depth];
  float tmaxlist[max_tree_depth];
  int axis;
  float dist0, dist1;

  if (!aabb_intersect(b->aabb, ray, invdir, &tmin, &tmax)) return 1;

  while (1)
  {
    while (1)
    {
      axis = node->data & 3;
      if(axis == 3) break;
      dist0 = (node->clip[below[axis]] - ray->pos[axis]) * invdir[axis];
      dist1 = (node->clip[far[axis]] - ray->pos[axis]) * invdir[axis];
      fbSide[0] = fbSide[1] = b->tree + (node->data >> 3);//bih_children(node, b);
      fbSide[far[axis]]++;

      if(node->data & 4)//bih_isClipNode(node))
      {
        tmin = max(tmin, dist0);
        tmax = min(tmax, dist1);
        if (tmin > tmax)
        {
          if ( nListPos > 0 )
          {
            --nListPos;
            tmin = tminlist[nListPos];
            tmax = min(tmaxlist[nListPos], 1.0f);
            node = nList[nListPos];
          }
          else return 1;
        }
        else node = b->tree + (node->data >> 3);//bih_children(node, b);
        continue;
      }

      if ( dist0 < tmin )
      {
        if ( dist1 > tmax )
        {
          if (nListPos > 0)
          {
            --nListPos;
            tmin = tminlist[nListPos];
            tmax = min(tmaxlist[nListPos], 1.0f);
            node = nList[nListPos];
          }
          else return 1;
        }
        else
        {
          tmin = max(tmin, dist1);
          node = fbSide[1];
        }
      }
      else if ( dist1 > tmax )
      {
        tmax = min(tmax, dist0);
        node = fbSide[0];
      }
      else
      {
        tminlist[nListPos] = max(tmin, dist1);
        tmaxlist[nListPos] = tmax;
        nList[nListPos] = fbSide[1];
        ++nListPos;

        node = fbSide[0];
        tmax = min(tmax, dist0);
      }
    }

    unsigned int tmp_index = node->data >> 3;//bih_prims(node);
    for (unsigned int i = 0; i < *(unsigned int*)node; i++)//bih_num(node); i++)
    {
      if(prims_intersect_visible(b->prims, b->prims->primid[tmp_index], ray)) return 0;
      tmp_index++;
    }

    while(1)
    {
      if (nListPos == 0) return 1;
      --nListPos;
      tmin = tminlist[nListPos];
      if (1.0f >= tmin)
      {
        tmax = min(tmaxlist[nListPos], 1.0f);
        node = nList[nListPos];
        break;
      }
    }
  }
}
