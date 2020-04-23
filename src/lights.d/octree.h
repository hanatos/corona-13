#pragma once

// XXX
// should we store a full cut? can we always simplify it if some nodes
// on the way are insignificant? what if only the leaf node is significant
// but the intermediate nodes not?
//
// advantages of full cut:
// can pull in new nodes during learning phase without storing a full list of
// new samples and their contributions
//
// XXX
// alternative:
// keep 5-10% for uniform lh
// store only important nodes in tree, not a complete cut
// discard least relevant ones without replacement
//
// XXX
// hybrid:
// store full cut but allow overlapping entry points in cdf:
// pdf is then sum of parent cut
// check if parent is in cut:
// - array is sorted:
// parent (or grand parent) would come directly before in array?
// - sum of all point pdfs (early out with aabb)

// node in adaptive octree that learns
// cuts through the light hierarchy qbvh
// per node
#define LH_CUTSIZE 32
typedef struct lh_onode_t
{
  float p [LH_CUTSIZE];    // accumulated power per node
  float p2[LH_CUTSIZE];    // p squared accumulated
  uint32_t n[LH_CUTSIZE];  // number of samples contributing to that
  // stores node indices to light hierarchy node buffer
  // sorted in order. -1 marks end of list.
  uint32_t cut[LH_CUTSIZE];
}
lh_onode_t;

static inline void
lh_init_cut(
    lh_onode_t *node,
    const lh_t *lh)
{
  // TODO: traverse light hierarchy depth first
  // and remember 32 (level 5) lights
}

// 
  // TODO: go through list, compute probability for given shading
  // point x ? directly sample?
//     const float *x,
//     float       *p)

static inline void
lh_calc_p(
    lh_onode_t  *node,
    const lh_t  *lh)
{
  // TODO: something that computes p from prior + accumulated stuff + variance?
}


// perform one iteration of cut refinement.  will split the most important node
// reference if the cut still has enough room. else it will take the most
// insignificant ones and replace by the common parent node
static inline void
lh_refine_cut(
    lh_onode_t *node,
    const lh_t *lh)
{
  // float avg_p = 0.0f;
  int cnt = 0;
  float pm =  FLT_MAX, pM = -FLT_MAX;
  uint32_t im = -1, iM = -1;
  for(int i=0;i<LH_CUTSIZE;i++)
  {
    if(node->cut[i] == -1) break;
    // avg_p += node->p[i];
    cnt++;
    if(node->p[i] < pm)
    {
      pm = node->p[i];
      im = i;
    }
    if(node->p[i] > pM)
    {
      pM = node->p[i];
      iM = i;
    }
  }
  // avg_p /= cnt;

  if(pM > 4.0 * pm)
  {
    // if enough room for split (+3):
    if(cnt < LH_CUTSIZE - 4)
    {
    // split: easy
    // remove iM from list and replace by its <= four children if it has any
    // init p and n by p/4 and n/4
    }
    else
    {
    // else: make some room:
    // merge: slightly harder:
    // look to the left and right of im and find nodes with same parent
    // collect all of these and replace by common parent
    // init p and n by sum(p) and sum(n)
    }
  }
}
