// lbvh cannibalised beyond recongnition from pbrt
// actually the "l" only stands for "lame" now.
#include <float.h>

#define LBVH_LEAF_PRIMS 7

typedef struct lbvh_node_t
{
  float aabb[2][6];  // two child boxes
  uint64_t child[2]; // -(prim<<5|num_prims) or child index
}
lbvh_node_t;

typedef struct lbvh_t
{
  float root_aabb[6]; // bounding box of whole thing
  float *aabb;        // tmp memory for build: bounding boxes
  lbvh_node_t *nodes;
  uint32_t *primid;   // primitive id permutation
  uint32_t num_prims;
  uint32_t max_nodes;
  uint32_t num_nodes;
  uint32_t num_alloced_prims;
}
lbvh_t;

static inline void lbvh_init(lbvh_t *b, int num_prims)
{
  // alloc morton and aabb and nodes and num_prims
  b->num_prims = num_prims;
  b->num_alloced_prims = num_prims;
  b->aabb = malloc(sizeof(float)*6*num_prims);
  b->primid = malloc(sizeof(*b->primid)*num_prims);
  b->max_nodes = MAX(2, 2*num_prims+1);
  b->num_nodes = 0;
  b->nodes = malloc(sizeof(lbvh_node_t)*b->max_nodes);
}

static inline void lbvh_cleanup(lbvh_t *b)
{
  free(b->aabb);
  free(b->nodes);
  free(b->primid);
  memset(b, 0, sizeof(*b));
}

static inline void lbvh_make_leaf(
    lbvh_t *b,
    int parent,
    int child,
    int start,
    int end)
{
  const int num = end-start;
  assert(num < (1<<5));
  assert(num > 0);
  assert(start < b->num_prims);
  assert(start >= 0);
  b->nodes[parent].child[child] = (1ul<<63) | ((start<<5)|num);
  for(int k=0;k<3;k++) b->nodes[parent].aabb[child][k] =  FLT_MAX;
  for(int k=3;k<6;k++) b->nodes[parent].aabb[child][k] = -FLT_MAX;
  for(int i=start;i<end;i++)
  { // expand bounding box
    const int ii = b->primid[i];
    assert(ii < b->num_prims);
    assert(ii >= 0);
    for(int k=0;k<3;k++) b->nodes[parent].aabb[child][k] =  MIN(
        b->nodes[parent].aabb[child][k], b->aabb[6*ii+k]);
    for(int k=3;k<6;k++) b->nodes[parent].aabb[child][k] =  MAX(
        b->nodes[parent].aabb[child][k], b->aabb[6*ii+k]);
  }
}

static inline void lbvh_build_median_rec(
    lbvh_t *b,
    int parent,
    int child,
    int start,
    int end)
{
  if(end-start < LBVH_LEAF_PRIMS) return lbvh_make_leaf(b, parent, child, start, end);
  int n = b->num_nodes++;
  assert(n < b->max_nodes);
  float aabb[6] = {FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};
  for(int i=start;i<end;i++)
  {
    for(int k=0;k<3;k++) aabb[k] = MIN(aabb[k], b->aabb[6*b->primid[i]+k]);
    for(int k=3;k<6;k++) aabb[k] = MAX(aabb[k], b->aabb[6*b->primid[i]+k]);
  }
  int dim = 0;
  if(aabb[4]-aabb[1] > aabb[dim+3]-aabb[dim]) dim = 1;
  if(aabb[5]-aabb[2] > aabb[dim+3]-aabb[dim]) dim = 2;
  float splitf = (aabb[dim+3]+aabb[dim])/2.0f;
  int split = end;
  for(int i=start;i<split;i++)
  {
    if((b->aabb[6*b->primid[i]+dim+3]+b->aabb[6*b->primid[i]+dim])/2.0f > splitf)
    {
      uint32_t tmp = b->primid[i];
      b->primid[i--] = b->primid[--split];
      b->primid[split] = tmp;
    }
  }
  if(end == split || split == start)
    split = (start+end)/2;
  assert(end > split);
  assert(split > start);
  lbvh_build_median_rec(b, n, 0, start, split);
  lbvh_build_median_rec(b, n, 1, split, end);
  if(parent >= 0)
  {
    for(int k=0;k<3;k++) b->nodes[parent].aabb[child][k] = MIN(
        b->nodes[n].aabb[0][k],
        b->nodes[n].aabb[1][k]);
    for(int k=3;k<6;k++) b->nodes[parent].aabb[child][k] = MAX(
        b->nodes[n].aabb[0][k],
        b->nodes[n].aabb[1][k]);
    // init child pointer
    b->nodes[parent].child[child] = n;
  }
}

static inline float _lbvh_dist(
    const float *aabb,
    const float *x)
{
  float dist = 0.0f;
  for(int k=0;k<3;k++)
  {
    float d = MAX(fabsf(x[k]-(aabb[3+k]+aabb[k])/2.0f) - (aabb[3+k]-aabb[k])/2.0f, 0.0f);
    dist += d*d;
  }
  return dist;
}

static inline void lbvh_build(lbvh_t *b)
{
  // compute root bounding box:
  for(int k=0;k<3;k++) b->root_aabb[k] =  FLT_MAX;
  for(int k=3;k<6;k++) b->root_aabb[k] = -FLT_MAX;
  for(int i=0;i<b->num_prims;i++)
  {
    for(int k=0;k<3;k++) b->root_aabb[k] = MIN(b->root_aabb[k], b->aabb[6*i+k]);
    for(int k=3;k<6;k++) b->root_aabb[k] = MAX(b->root_aabb[k], b->aabb[6*i+k]);
  }
  for(int i=0;i<b->num_prims;i++)
    b->primid[i] = i;

  if(b->num_prims < 32)
  {
    b->num_nodes = 1;
    // fprintf(stderr, "[lbvh] making root a leaf!\n");
    b->nodes[0].child[0] = (1ul<<63) | b->num_prims;
    b->nodes[0].child[1] = (1ul<<63);
    for(int c=0;c<2;c++)
    {
      for(int k=0;k<3;k++) b->nodes[0].aabb[c][k]   =  FLT_MAX;
      for(int k=0;k<3;k++) b->nodes[0].aabb[c][k+3] = -FLT_MAX;
    }
    for(int i=0;i<b->num_prims;i++)
    { // expand bounding box
      const int ii = b->primid[i];
      assert(ii < b->num_prims);
      assert(ii >= 0);
      for(int k=0;k<3;k++) b->nodes[0].aabb[0][k] =  MIN(
                           b->nodes[0].aabb[0][k], b->aabb[6*ii+k]);
      for(int k=3;k<6;k++) b->nodes[0].aabb[0][k] =  MAX(
                           b->nodes[0].aabb[0][k], b->aabb[6*ii+k]);
    }
    return;
  }
  //const int bit_index = 60;
  b->num_nodes = 0;
  lbvh_build_median_rec(b, -1, -1, 0, b->num_prims);
}

// refit bvh based on updated leaf aabbs
static inline uint32_t lbvh_refit(lbvh_t *b, lbvh_node_t *node)
{
  // ideally this would catch all `root node is a leaf' cases:
  if(b->num_prims == 0) return 0;
  uint32_t num_prims = 0;

  // update root bounding box:
  for(int k=0;k<3;k++) b->root_aabb[k] =  FLT_MAX;
  for(int k=3;k<6;k++) b->root_aabb[k] = -FLT_MAX;
  for(int i=0;i<b->num_prims;i++)
  {
    for(int k=0;k<3;k++) b->root_aabb[k] = MIN(b->root_aabb[k], b->aabb[6*i+k]);
    for(int k=3;k<6;k++) b->root_aabb[k] = MAX(b->root_aabb[k], b->aabb[6*i+k]);
  }
  for(int c=0;c<2;c++)
  {
    if(node->child[c] & (1ul<<63))
    { // leaf
      for(int k=0;k<3;k++) node->aabb[c][k] =  FLT_MAX;
      for(int k=3;k<6;k++) node->aabb[c][k] = -FLT_MAX;
      const int prim =((1ul<<63) ^ node->child[c])>>5;
      const int cnt = node->child[c] & ~-(1<<5);
      num_prims += cnt;
      assert(prim < b->num_prims);
      assert(prim + cnt <= b->num_prims);
      for(int i=prim;i<prim+cnt;i++)
      { // expand bounding box
        const int ii = b->primid[i];
        for(int k=0;k<3;k++) node->aabb[c][k] =  MIN(
                             node->aabb[c][k], b->aabb[6*ii+k]);
        for(int k=3;k<6;k++) node->aabb[c][k] =  MAX(
                             node->aabb[c][k], b->aabb[6*ii+k]);
      }
    }
    else
    { // inner
      assert(node->child[c] < b->num_nodes);
      assert(node->child[c] < b->max_nodes);
      lbvh_node_t *child = b->nodes + node->child[c];
      num_prims += lbvh_refit(b, child);
      for(int k=0;k<3;k++) node->aabb[c][k] = MIN(child->aabb[0][k], child->aabb[1][k]);
      for(int k=3;k<6;k++) node->aabb[c][k] = MAX(child->aabb[0][k], child->aabb[1][k]);
    }
  }
  if(node == b->nodes) assert(num_prims == b->num_prims);
  return num_prims;
}

#undef LBVH_LEAF_PRIMS
