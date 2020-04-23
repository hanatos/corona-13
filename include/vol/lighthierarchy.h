#pragma once

#include "vol/vol.h"
#include "vol/shaders.h"
#include "vol/interpolation.h"
#include "spectrum.h"

static inline int _vol_sample_binary(const float p, float *rand, float *pdf)
{
  assert(p == p);
  assert(*rand == *rand);
  if(p < 1 && *rand >= p)
  {
    *rand = (*rand-p)/(1.0f-p);
    *pdf *= 1.0f-p;
    assert(*pdf > 0.0);
    return 1;
  }
  *rand /= p;
  *pdf  *= p;
  assert(*pdf > 0.0);
  return 0;
}

static inline void _eval_culln(
    const float *cullx, // in object space
    const float *culln, // in object space
    const int x,
    const int y,
    const int z,
    const float voxel_size,
    const int l, // level in tree, 0..3*(tree depth)-1
    float *v)
{
  if(l<3 || !culln || !cullx)
  {
    for(int k=0;k<8;k++) v[k] = 1.0f;
    return;
  }
  // compute 8 cosines <(p - cullx),culln>/|p-cullx|
  const float vs = voxel_size * (1<<l);
  for(int i=0;i<8;i++)
  {
    float p[3] = {
      x+.25f + ((i&1)?.5f:0.f),
      y+.25f + ((i&2)?.5f:0.f),
      z+.25f + ((i&4)?.5f:0.f)
    };
    for(int k=0;k<3;k++) p[k] *= vs;
    const float d[3] = {
      p[0]-cullx[0], p[1]-cullx[1], p[2]-cullx[2] };
    // this would really be 1 + cos, but maybe we can hard-cut off a bit:
    v[i] = MAX(0.0, 1.05f + dotproduct(d, culln)/sqrtf(dotproduct(d, d)));
  }
}

static inline void _eval_cullp(
    const float *cullv,
    const uint16_t *binp,
    float *binary)
{
  // 0 | 1 2 | 3 4 5 6
  float v[8]; // convert binary search to voxel probabilities:
  v[0] = binp[0]/(float)0x10000 * binp[1]/(float)0x10000 * binp[3]/(float)0x10000;
  v[1] = binp[0]/(float)0x10000 * binp[1]/(float)0x10000 * (1.0f-binp[3]/(float)0x10000);
  v[2] = binp[0]/(float)0x10000 * (1.0f-binp[1]/(float)0x10000) * binp[4]/(float)0x10000;
  v[3] = binp[0]/(float)0x10000 * (1.0f-binp[1]/(float)0x10000) * (1.0f-binp[4]/(float)0x10000);
  v[4] = (1.0f-binp[0]/(float)0x10000) * binp[2]/(float)0x10000 * binp[5]/(float)0x10000;
  v[5] = (1.0f-binp[0]/(float)0x10000) * binp[2]/(float)0x10000 * (1.0f-binp[5]/(float)0x10000);
  v[6] = (1.0f-binp[0]/(float)0x10000) * (1.0f-binp[2]/(float)0x10000) * binp[6]/(float)0x10000;
  v[7] = (1.0f-binp[0]/(float)0x10000) * (1.0f-binp[2]/(float)0x10000) * (1.0f-binp[6]/(float)0x10000);
  // combine with culling values:
  for(int k=0;k<8;k++) v[k] *= cullv[k];
  // convert back to binary search tree:
  float total = v[0]+v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7], t = 0.;
  memset(binary, 0, sizeof(float)*7);
  if(total > 0.0)
  {
    binary[0] = (v[0]+v[1]+v[2]+v[3])/total; // p(z < .5):
    t = v[0]+v[1]+v[2]+v[3];
    if(t > 0.0) binary[1] = (v[0] + v[1])/t; // p(y < .5 | z  < .5)
    t = v[4]+v[5]+v[6]+v[7];
    if(t > 0.0) binary[2] = (v[4] + v[5])/t; // p(y < .5 | z >= .5)

    t = v[0]+v[1];
    if(t > 0.0) binary[3] = v[0]/t; // p(x < .5 | y <  .5, z <  .5)
    t = v[2]+v[3];
    if(t > 0.0) binary[4] = v[2]/t; // ..         y >= .5, z <  .5
    t = v[4]+v[5];
    if(t > 0.0) binary[5] = v[4]/t; //            y <  .5, z >= .5
    t = v[6]+v[7];
    if(t > 0.0) binary[6] = v[6]/t; //            y >= .5, z >= .5
  }
}

// traverse light hierarchy, fill position, return L_e (shader evaluation in target voxel * density)
static inline mf_t vol_lighthierarchy_sample_point(
    const vol_tree_t *const tree,
    const mf_t lambda,
    const float time,
    const int lod, // level of detail, 0 means finest, 1 means 8x coarser etc
    const int use_max, 
    float rand0,
    float rand1,
    float rand2,
    const float *cullx_ws,
    const float *culln_ws,
    float *pos,
    mf_t  *pdf_out)
{
  const vol_node_t *n = tree->root_node;
  const vol_payload_t *payload = tree->root_payload;
  const vol_payload_flux_t *node = tree->root_light;
  int isstatic = tree->header->isstatic;

  const int max_depth = tree->header->depth;
  int depth = 0;

  float cullxs[3], *cullx = 0;
  float cullns[3], *culln = 0;
  if(cullx_ws && culln_ws)
  {
    for(int k=0;k<3;k++) cullxs[k] = cullx_ws[k];
    for(int k=0;k<3;k++) cullns[k] = culln_ws[k];
    vol_transform_o2w(tree, cullxs, 1); // transform object to world space
    vol_transform_o2w(tree, cullns, 0); // det == 1, so transforming normals like that is okay
    for(int k=0;k<3;k++)
      cullxs[k] -= tree->aabb[k];
    cullx = cullxs;
    culln = cullns;
  }

  // remember pdf of what we've done
  float pdf = 1.0f;
  // track center point
  int x=0, y=0, z=0;
  float L_e = 0.0f;
  float cullv[8];
  float binary[7];

  // fprintf(stderr, "sample\n");
  while(1)
  {
    // sample one child a from base block node->p[8][.]
    int k = 0, j = 0, i = 0;
    const uint16_t (*data)[7] = use_max ? node->q : node->p;
    // fprintf(stderr, "level %d x %u %u %u\n", 3*max_depth-3*depth-0, x, y, z);
    _eval_culln(cullx, culln, x, y, z, tree->voxel_size, 3*max_depth-3*depth-0, cullv);
    _eval_cullp(cullv, data[8], binary);
    k = _vol_sample_binary(binary[0      ], &rand2, &pdf);
    j = _vol_sample_binary(binary[1+k    ], &rand1, &pdf);
    i = _vol_sample_binary(binary[3+2*k+j], &rand0, &pdf);
    // k = _vol_sample_binary(data[8][0      ]/(float)0x10000, &rand2, &pdf);
    // j = _vol_sample_binary(data[8][1+k    ]/(float)0x10000, &rand1, &pdf);
    // i = _vol_sample_binary(data[8][3+2*k+j]/(float)0x10000, &rand0, &pdf);

    const int c = (k<<2)|(j<<1)|i; // child index
    // child index in full node
    int child = (k<<8)|(j<<5)|(i<<2);
    x = (x<<1)|i; y = (y<<1)|j; z = (z<<1)|k;

    // redo for child voxel:
    _eval_culln(cullx, culln, x, y, z, tree->voxel_size, 3*max_depth-3*depth-1, cullv);
    _eval_cullp(cullv, data[c], binary);
    k = _vol_sample_binary(binary[0      ], &rand2, &pdf);
    j = _vol_sample_binary(binary[1+k    ], &rand1, &pdf);
    i = _vol_sample_binary(binary[3+2*k+j], &rand0, &pdf);
    // k = _vol_sample_binary(data[c][0      ]/(float)0x10000, &rand2, &pdf);
    // j = _vol_sample_binary(data[c][1+k    ]/(float)0x10000, &rand1, &pdf);
    // i = _vol_sample_binary(data[c][3+2*k+j]/(float)0x10000, &rand0, &pdf);

    child |= (k<<7)|(j<<4)|(i<<1);
    x = (x<<1)|i; y = (y<<1)|j; z = (z<<1)|k;

    // redo for actual voxel data obtained from regular tree:
    // get flux of eight children in full payload in tree, normalise, sample third level
    float L[8] = {0.0f};
    _eval_culln(cullx, culln, x, y, z, tree->voxel_size, 3*max_depth-3*depth-2, cullv);
    for(int ii=0;ii<8;ii++)
    {
      const int idx = child + (((ii&4)<<4)|((ii&2)<<2)|(ii&1));
      float result[2] = {0.0f};
      vol_payload_get(payload, idx, time, isstatic, result);
      if(depth+1 < max_depth)
        L[ii] = cullv[ii] * (use_max ? result[0] : result[1]);
      else if(result[0] > 0.0f && result[1] > 0.0f)
        L[ii] = cullv[ii] * result[0] * mf(tree->shader(result[0], result[1], lambda),0);
    }
    const float total = L[0]+L[1]+L[2]+L[3]+L[4]+L[5]+L[6]+L[7];
    if(total <= 0.0) return mf_set1(0.0f);
    k = _vol_sample_binary((L[0]+L[1]+L[2]+L[3])/total, &rand2, &pdf);
    j = _vol_sample_binary((L[4*k]+L[4*k+1])/(L[4*k]+L[4*k+1]+L[4*k+2]+L[4*k+3]), &rand1, &pdf);
    i = _vol_sample_binary(L[4*k+2*j]/(L[4*k+2*j]+L[4*k+2*j+1]), &rand0, &pdf); 
    child |= (k<<6)|(j<<3)|i;
    x = (x<<1)|i; y = (y<<1)|j; z = (z<<1)|k;
    L_e = L[k*4+j*2+i];

    if(++depth >= max_depth-lod) break;
    payload = vol_node_get_payload(tree, n, child);
    if(!payload) return mf_set1(0.0f); // sampled empty location, bummer. this should only happen during debugging!
    assert(payload);
    isstatic = vol_node_child_static(n, child);
    node = vol_node_get_flux(tree, n, child);
    n = vol_node_get_child(tree, n, child);
  }
  // sample uniformly inside last voxel, adjust pdf
  // clamp random numbers to avoid rounding issues when asking for the pdf:
  // if we overflow to pos[0]+1 for instance, the returned pdf will be different.
  pos[0] = x + CLAMP(rand0, 0, 0.999);
  pos[1] = y + CLAMP(rand1, 0, 0.999);
  pos[2] = z + CLAMP(rand2, 0, 0.999);
  // adjust voxel size for max depth/lod
  const float vs = tree->voxel_size * (1<<(3*lod));
  pdf *= 1.0f/(vs*vs*vs);

  // transform grid to object space:
  for(int k=0;k<3;k++) pos[k] = pos[k]*vs + tree->aabb[k];

  for(int k=0;k<3;k++)
    assert(pos[k] >= tree->aabb[k] && pos[k] <= tree->aabb[k+3]);

  vol_transform_o2w(tree, pos, 1); // transform object to world space

  if(pdf_out) *pdf_out = mf_set1(pdf);
  assert(pdf > 0.0f);
  assert(L_e > 0.0f);
  return mf_set1(L_e);
}

// pdf of light hierarchy point sampling.
static inline mf_t vol_lighthierarchy_pdf_point(
    const vol_tree_t *const tree,
    const mf_t lambda,
    const float time,
    const int lod, // level of detail, 0 means finest, 1 means 8x coarser etc
    const int use_max,
    const float *cullx_ws,
    const float *culln_ws,
    const float *pos_ws)
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  vol_transform_w2o(tree, pos, 1); // transform to object space
  if((pos[0] < tree->content_box[0] || pos[0] > tree->content_box[3]) ||
     (pos[1] < tree->content_box[1] || pos[1] > tree->content_box[4]) ||
     (pos[2] < tree->content_box[2] || pos[2] > tree->content_box[5]))
  {
    return mf_set1(0.f);
  }
  const int x = (pos[0] - tree->aabb[0])/tree->voxel_size; // transform to voxel index space
  const int y = (pos[1] - tree->aabb[1])/tree->voxel_size;
  const int z = (pos[2] - tree->aabb[2])/tree->voxel_size;
  // now in grid coordinates

  const vol_node_t *n = tree->root_node;
  const vol_payload_t *payload = tree->root_payload;
  const vol_payload_flux_t *node = tree->root_light;
  int isstatic = tree->header->isstatic;

  const int max_depth = tree->header->depth;
  int depth = 0;

  float pdf = 1.0f;
  int l = max_depth*3-1;
  int k = 0, j = 0, i = 0;

  float cullv[8];
  float binary[7];
  float cullxs[3], *cullx = 0;
  float cullns[3], *culln = 0;
  if(cullx_ws && culln_ws)
  {
    for(int k=0;k<3;k++) cullxs[k] = cullx_ws[k];
    for(int k=0;k<3;k++) cullns[k] = culln_ws[k];
    vol_transform_o2w(tree, cullxs, 1); // transform object to world space
    vol_transform_o2w(tree, cullns, 0); // det == 1, so transforming normals like that is okay
    for(int k=0;k<3;k++)
      cullxs[k] -= tree->aabb[k];
    cullx = cullxs;
    culln = cullns;
  }

  // fprintf(stderr, "pdf\n");
  while(l>=0)
  {
    const uint16_t (*data)[7] = use_max ? node->q : node->p;
    i = (x>>l)&1; j = (y>>l)&1; k = (z>>l)&1; l--;
    // fprintf(stderr, "level %d x %u %u %u\n", l+2, x>>(l+2), y>>(l+2), z>>(l+2));
    _eval_culln(cullx, culln, x>>(l+2), y>>(l+2), z>>(l+2), tree->voxel_size, l+2, cullv);
    _eval_cullp(cullv, data[8], binary);
    if(k) pdf *= 1.0f-binary[0      ];//data[8][0      ]/(float)0x10000;
    else  pdf *=      binary[0      ];//data[8][0      ]/(float)0x10000;
    if(j) pdf *= 1.0f-binary[1+k    ];//data[8][1+k    ]/(float)0x10000;
    else  pdf *=      binary[1+k    ];//data[8][1+k    ]/(float)0x10000;
    if(i) pdf *= 1.0f-binary[3+2*k+j];//data[8][3+2*k+j]/(float)0x10000;
    else  pdf *=      binary[3+2*k+j];//data[8][3+2*k+j]/(float)0x10000;
    const int c = (k<<2)|(j<<1)|i;
    int child = (k<<8)|(j<<5)|(i<<2);

    i = (x>>l)&1; j = (y>>l)&1; k = (z>>l)&1; l--;
    _eval_culln(cullx, culln, x>>(l+2), y>>(l+2), z>>(l+2), tree->voxel_size, l+2, cullv);
    _eval_cullp(cullv, data[c], binary);
    if(k) pdf *= 1.0f-binary[0      ];//data[c][0      ]/(float)0x10000;
    else  pdf *=      binary[0      ];//data[c][0      ]/(float)0x10000;
    if(j) pdf *= 1.0f-binary[1+k    ];//data[c][1+k    ]/(float)0x10000;
    else  pdf *=      binary[1+k    ];//data[c][1+k    ]/(float)0x10000;
    if(i) pdf *= 1.0f-binary[3+2*k+j];//data[c][3+2*k+j]/(float)0x10000;
    else  pdf *=      binary[3+2*k+j];//data[c][3+2*k+j]/(float)0x10000;
    child |= (k<<7)|(j<<4)|(i<<1);
    if (pdf <= 0.0f) return mf_set1(0.0f);

    i = (x>>l)&1; j = (y>>l)&1; k = (z>>l)&1; l--;
    float L[8] = {0.0f};
    _eval_culln(cullx, culln, x>>(l+2), y>>(l+2), z>>(l+2), tree->voxel_size, l+2, cullv);
    for(int ii=0;ii<8;ii++)
    {
      const int idx = child + (((ii&4)<<4)|((ii&2)<<2)|(ii&1));
      float result[2] = {0.0f};
      vol_payload_get(payload, idx, time, isstatic, result);
      if(depth+1 < max_depth)
        L[ii] = cullv[ii]*(use_max ? result[0] : result[1]);
      else if(result[0] > 0.0f && result[1] > 0.0f)
        // XXX FIXME: we actually need to propagate everything to all wavelengths here!
        L[ii] = cullv[ii] * result[0] * mf(tree->shader(result[0], result[1], lambda),0);
    }
    const float total = L[0]+L[1]+L[2]+L[3]+L[4]+L[5]+L[6]+L[7];
    const float sub_totalj = (L[4*k]+L[4*k+1]+L[4*k+2]+L[4*k+3]);
    const float sub_totali = (L[4*k+2*j]+L[4*k+2*j+1]);
    if(total <= 0.0 || sub_totalj <= 0.0 || sub_totali <= 0.0) return mf_set1(0.0f);
    if(k) pdf *= 1.0f-(L[0]+L[1]+L[2]+L[3])/total;
    else  pdf *=      (L[0]+L[1]+L[2]+L[3])/total;
    if(j) pdf *= 1.0f-(L[4*k]+L[4*k+1])/sub_totalj;
    else  pdf *=      (L[4*k]+L[4*k+1])/sub_totalj;
    if(i) pdf *= 1.0f-L[4*k+2*j]/sub_totali;
    else  pdf *=      L[4*k+2*j]/sub_totali;
    child |= (k<<6)|(j<<3)|i;
    assert(pdf == pdf);

    if(++depth >= max_depth-lod) break;
    payload = vol_node_get_payload(tree, n, child);
    if(!payload) return mf_set1(0.0f); // sampled empty location, bummer. this should only happen during debugging!
    assert(payload);
    isstatic = vol_node_child_static(n, child);
    node = vol_node_get_flux(tree, n, child);
    n = vol_node_get_child(tree, n, child);
  }
  // adjust for LOD/max level
  const float vs = tree->voxel_size * (1<<(3*lod));
  pdf *= 1.0f/(vs*vs*vs);
  assert(pdf == pdf);
  return mf_set1(pdf);
}

#if 0
#include "vol/trace_octree.h"
// trace through light hierarchy and find out pdf of sampling the given direction
static inline float vol_lighthierarchy_pdf_direction(
    const vol_tree_t *tree,              // volume tree
    const float *pos_ws,                 // world space position of query ray
    const float *dir_ws,                 // world space direction of query ray
    const float lambda,                  // wavelength
    float time,                          // time for motion blur access
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    const int lod2,                       // LOD: 0 is finest, 1 is 8x coarser, ..
    const int use_max,
    const float *cullx_ws,
    const float *culln_ws)
{
  const int lod = 0;
  const float max_dist = FLT_MAX;
  float accum_pdf = 0.0f;
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  float cullxs[3], *cullx = 0;
  float cullns[3], *culln = 0;
  if(cullx_ws && culln_ws)
  {
    for(int k=0;k<3;k++) cullxs[k] = cullx_ws[k];
    for(int k=0;k<3;k++) cullns[k] = culln_ws[k];
    vol_transform_o2w(tree, cullxs, 1); // transform object to world space
    vol_transform_o2w(tree, cullns, 0); // det == 1, so transforming normals like that is okay
    for(int k=0;k<3;k++)
      cullxs[k] -= tree->aabb[k];
    cullx = cullxs;
    culln = cullns;
  }
#define VOL_TRACE_INIT
#define VOL_TRACE_LEAF \
  const float tc = .5f*(tmax + tmin);\
  if(tmax > tmin)\
  {\
    const float posv[3] = { pos_ws[0] + tc*dir_ws[0], pos_ws[1] + tc*dir_ws[1], pos_ws[2] + tc*dir_ws[2] };\
    const float pdf = vol_lighthierarchy_pdf_point(tree, lambda, time, lod, 1, cullx_ws, culln_ws, posv);\
    accum_pdf += (tmax*tmax*tmax - tmin*tmin*tmin)/3.0f * pdf;\
  }
#define VOL_TRACE_MISS return 0.0f;
#define VOL_TRACE_RETURN \
  const float accum_pdf2 = vol_lighthierarchy_pdf_direction2(tree, pos_ws, dir_ws, lambda, time, interpolation, lod, use_max, cullx_ws, culln_ws);\
  fprintf(stderr, "pdfs: %g %g\n", accum_pdf2, accum_pdf);\
  return accum_pdf;
#include "vol/trace_impl.inc"
}
#else
#include "vol/trace_octree.h"
#endif

// TODO: remove this with the sample point func above
static inline float vol_lighthierarchy_sample_point_fnee(
    const vol_tree_t *const tree,
    const float lambda,
    const float time,
    const int lod, // level of detail, 0 means finest, 1 means 8x coarser etc
    float rand0,
    float rand1,
    float rand2,
    float *pos,
    float *pdf_out)
{
  const float *b = tree->content_box;
  pos[0] = b[0] + rand0 * (b[3]-b[0]);
  pos[1] = b[1] + rand1 * (b[4]-b[1]);
  pos[2] = b[2] + rand2 * (b[5]-b[2]);
  vol_transform_o2w(tree, pos, 1); // transform object to world space
  return 0.f;
}

static inline float vol_lighthierarchy_pdf_fnee_direction(
    const vol_tree_t *const tree,
    const float lambda,
    float time,
    const int lod,                       // level of detail, 0 means finest, 1 means 8x coarser etc
    const float *pos_ws,
    const float *omega_ws)
{
  float omega[3] = {omega_ws[0], omega_ws[1], omega_ws[2]};
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  vol_transform_w2o(tree, omega, 0);
  vol_transform_w2o(tree, pos, 1);

  // intersect with content_box
  const float *cbox = tree->content_box;
  float invdir[3];
  for(int k=0;k<3;k++) invdir[k] = 1.0f/omega[k];

  // intersect bounding box, get entry and exit distances tmin, tmax
  float tmin = 0.0f;
  float tmax = FLT_MAX;
  for(int k=0;k<3;k++)
  {
    const float t0 = (cbox[k+3] - pos[k])*invdir[k];
    const float t1 = (cbox[k]   - pos[k])*invdir[k];
    if(t0 <= t1)
    {
      tmin = t0 > tmin ? t0 : tmin;
      tmax = t1 < tmax ? t1 : tmax;
    }
    else
    {
      tmin = t1 > tmin ? t1 : tmin;
      tmax = t0 < tmax ? t0 : tmax;
    }
  }
  // miss bounding box?
  if(tmin > tmax)
  {
    return 0.f;
  }

  const float volume = (tree->content_box[3]-tree->content_box[0]) *
      (tree->content_box[4]-tree->content_box[1]) *
      (tree->content_box[5]-tree->content_box[2]);
  return (tmax*tmax*tmax-tmin*tmin*tmin) / (3.f*volume);
}
