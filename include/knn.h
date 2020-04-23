#ifndef _KNN_H
#define _KNN_H

#include "heap.h"

// find k nearest neighbours in a kd tree

#define KNN_DIM 5

typedef struct knn_point_t
{
  float pos[5];         // usually [x,y,z,lambda,time]
  uint32_t axis : 3;    // axis 0..KNN_DIM-1
  uint32_t right : 29;
  uint32_t payload;
}
knn_point_t;

static inline void _knn_aabb(float* aabb, const knn_point_t *points, int si, int ei)
{
  for(int k=0;k<KNN_DIM;k++) aabb[k+KNN_DIM] = aabb[k] = points[si].pos[k];
  for (int i = si+1; i < ei; i++)
  {
    for(int k=0;k<KNN_DIM;k++) aabb[k]         = fminf(points[i].pos[k], aabb[k]);
    for(int k=0;k<KNN_DIM;k++) aabb[k+KNN_DIM] = fmaxf(points[i].pos[k], aabb[k+KNN_DIM]);
  }
}

#define KNN_SWAP(a, b) {knn_point_t temp = points[a]; points[a] = points[b]; points[b] = temp;}

static inline int _knn_split_sort(knn_point_t *points, int si, int ei, int axis, float d)
{
  int m = 0; 
  while (si < ei)
  {
    while ((si < ei) && (points[si].pos[axis] < d)) si++;
    while ((si < ei) && (points[ei-1].pos[axis] >= d)) ei--;
    if (si < ei-1) KNN_SWAP(si, ei-1);
  }
  m = si;
  return m;
}

static inline void knn_build_rec(knn_point_t* points, int si, int ei)
{
  if(si == ei) return; // no samples here.
  if(si == ei - 1)
  { // only one sample left
    points[si].right = 0;
    points[si].axis = KNN_DIM;
    return;
  }

  // find longest axis and split it in the middle
  float bbox[2*KNN_DIM];
  _knn_aabb(bbox, points, si, ei);
  float temp[KNN_DIM]; for(int k=0;k<KNN_DIM;k++) temp[k] = bbox[k+KNN_DIM] - bbox[k];
  int axis = 0;
  for(int k=0;k<KNN_DIM;k++) if(temp[k] > temp[axis]) axis = k;
  float d = (bbox[axis] + bbox[KNN_DIM+axis])*0.5f;

  // find the sample closest to the split
  int pIdx = si;
  float dist = fabsf(d - points[si].pos[axis]);
  for (int i = si+1; i < ei; i++)
  {
    float npdist = fabsf(d - points[i].pos[axis]);
    if (npdist < dist)
    {
      dist = npdist;
      pIdx = i;
    }
  }

  d = points[pIdx].pos[axis]; // move the plane to the sample point
  KNN_SWAP(si, pIdx); // now swap the sample to the beginning of the list
  int m = _knn_split_sort(points, si+1, ei, axis, d);
  points[si].right = m;
  points[si].axis = axis;

  // TODO: benchmark when this actually makes sense to run in parallel
// #pragma omp task
  if (si+1 < m) knn_build_rec(points, si+1, m); // sort the left part
// #pragma omp task
  if (m < ei)   knn_build_rec(points, m, ei);   // sort the right part
  if(m >= ei) points[si].right = 0; // the right part is 0
}
#undef KNN_SWAP

// compute squared distance
static inline float _knn_dist(const float *x, const float *y)
{
  float dist = (x[0] - y[0])*(x[0] - y[0]);
  for(int k=1;k<KNN_DIM;k++)
    dist += (x[k] - y[k])*(x[k] - y[k]);
  return dist;
}

static inline int knn_find(
    knn_point_t *points,  // list of sample points in n-dimensions
    int num_points,       // number of points
    const float *q,       // query position in n-dimensions
    heap_t *res,          // result set stored in this heap. we'll find as many results as the heap can hold, or num_points.
    float mdist)          // maximum distance from query point.
{
  if(num_points == 0) return 1;
  int stack[128];
  int top = 0;
  stack[top] = 0;

  // uint64_t num_inserts = 0;
  // uint64_t num_steps = 0;
  // uint64_t num_unpruned = 0;
  // float mean_dist = 0;
  // clean up old results
  heap_clear(res);

  while(top >= 0)
  {
    const knn_point_t *p = points + stack[top];
    float dist = _knn_dist(p->pos, q);

    // num_steps ++;

    if(dist < mdist)
    {
      // num_inserts ++;
      uint64_t key;
      float val;
      if(heap_full(res)) heap_remove(res, &key, &val);
      heap_insert(res, stack[top], dist);
      if(heap_full(res))
      { // enable pruning only when we already found all k neighbours
        mdist = res->vals[0]; // head element is largest distance
        // mean_dist += mdist;
      }
      // else num_unpruned ++;
    }

    if(p->axis == KNN_DIM)
    { // check for leaf
      top--;
    }
    else
    {
      int ax = p->axis;
      float d = q[ax] - p->pos[ax];
      if(d < 0)
      { // left side first
        int left = stack[top] + 1;
        if ((p->right != 0) && (left != p->right) && (d*d < mdist))
          stack[top++] = p->right;
        stack[top] = left;
      }
      else
      { // right side first
        int left = stack[top] + 1;
        if ((left != p->right) && (d*d < mdist))
          stack[top++] = left;
        if (p->right != 0)
          stack[top] = p->right;
        else
          top--;
      }
    }
  }
  // fprintf(stderr, "[knn] steps %lu/%u inserts %lu mean dist %f unpruned %lu\n", num_steps, num_points, num_inserts, mean_dist/(num_inserts-num_unpruned), num_unpruned);
  // for(int k=0;k<res->end;k++)
  //   fprintf(stderr, "[knn] heap %d/%d %f\n", k, res->end, res->vals[k]);
  return 0;
}
#undef KNN_DIM


#endif
