#pragma once

#include <string.h>

#define LOW_BITS

// 5120 bytes
typedef struct vol_payload_compressed_t
{
  uint16_t d[512]; // master slice density
  uint16_t t[512]; // master slice temperature
#ifdef LOW_BITS
  uint8_t ref[VOL_MOTION_SAMPLES][16][3]; // offsets (6-bit each)
#else // pad up to 8-bits per ref for mental sanity and readout speed:
  uint8_t ref[VOL_MOTION_SAMPLES][64]; // offsets (8-bit each)
#endif
}
vol_payload_compressed_t;

typedef vol_payload_compressed_t vol_payload_t;

static inline size_t vol_payload_static_size()
{
  return sizeof(uint16_t)*1024;
}

// implement some sophisticated vector quantisation scheme to compress 4D grids.
// this is based on linde/buzo/gray and the fast variant
// "a fast mean-distance-ordered partial codebook search algorithm for image
//  vector quantization" (ra and kim 1993)

// write subsampled block reference index to compressed struct
static inline void _vol_payload_compress_entry(
    vol_payload_compressed_t *p,
    int i64, // subsampled block index (spatial) 0<=i64<64
    int t,   // time index 0<=t<VOL_MOTION_SAMPLES
    int r)   // desired reference block 0<r<64 in the master slice. =0 means empty block.
{
#ifdef LOW_BITS
  int k = i64 * 6;         // bit address of 6-bit reference
  int kb = k / 24;         // block index (one block has 3 bytes)
  int ki  = (k - 24*kb)/6; // bit address within one block / 6 bit per number
  switch(ki)
  {
    case 0: // [0]:6
      p->ref[t][kb][0] &= 0xc0;
      p->ref[t][kb][0] |= r;
      break;
    case 1: // [0]:2 and [1]:4
      p->ref[t][kb][0] &= 0x3f;
      p->ref[t][kb][0] |= (r&3)<<6;
      p->ref[t][kb][1] &= 0xf0;
      p->ref[t][kb][1] |= r>>2;
      break;
    case 2: // [1]:4 and [2]:2
      p->ref[t][kb][1] &= 0xf;
      p->ref[t][kb][1] |= (r&0xf)<<4;
      p->ref[t][kb][2] &= 0xfc;
      p->ref[t][kb][2] |= r>>4;
      break;
    default: // case 3: [2]:6
      p->ref[t][kb][2] &= 0x3;
      p->ref[t][kb][2] |= r<<2;
      break;
  }
#else
  p->ref[t][i64] = r;
#endif
}

// returns whether in fact the whole thing was static.
static inline int vol_payload_compress(
    const vol_payload_uncompressed_t *u,  // uncompressed grid
    vol_payload_t *p,                     // fill data into here
    int isstatic)                         // pass 1 when forcing to static
{
  if(isstatic)
  { // copy first time slice
    for(int k=0;k<512;k++) p->d[k] = float_to_half(CLAMP(u->d[k][0], 0, 65504.0));
    for(int k=0;k<512;k++) p->t[k] = float_to_half(CLAMP(u->t[k][0], 0, 65504.0));
    return 1;
  }
  const int num_it = 10;
  typedef struct codeword_t
  { // 2x2x2 block of our volume, for one time slice
    float d[8];
    float t[8];
  }
  __attribute__((aligned(32))) codeword_t;
  codeword_t codebook_old[64];
  codeword_t codebook_new[64];
  codeword_t uncompressed[VOL_MOTION_SAMPLES][64];
  // find initial codebook by selecting time slice with highest fill ratio:
  int master = 0, cnt = 0;
  int cnt_new[64];
  // ingest uncompressed block to simd-happy layout:
  for(int t=0;t<VOL_MOTION_SAMPLES;t++)
  {
    int new_cnt = 0;
    for(int i=0;i<512;i++)
    {
      const int i64 = ((i&0x6)>>1)|((i&0x30)>>2)|((i&0x180)>>3); // address of desired 4x4x4 block
      const int i8 = ((i&0x40)>>4)|((i&8)>>2)|(i&1);             // offset in 2x2x2 block
      uncompressed[t][i64].d[i8] = u->d[i][t];
      uncompressed[t][i64].t[i8] = u->t[i][t];
      if(u->d[i][t] > 0.0) new_cnt++;
    }
    if(new_cnt > cnt)
    {
      cnt = new_cnt;
      master = t;
    }
  }

  // init codebook to master time slice:
  memcpy(codebook_old, uncompressed[master], sizeof(codeword_t)*64);
  // and make codeword[0] empty to avoid special case
  memset(codebook_old, 0, sizeof(codeword_t));
  // default to empty blocks:
  for(int t=0;t<VOL_MOTION_SAMPLES;t++) for(int i=0;i<64;i++) _vol_payload_compress_entry(p, i, t, 0);

  // precompute half of squared mean distance:
  float mean[VOL_MOTION_SAMPLES][64] = {{0}};
  for(int t=0;t<VOL_MOTION_SAMPLES;t++) for(int i=0;i<64;i++)
    for(int j=0;j<8;j++)
      mean[t][i] += uncompressed[t][i].d[j] + 
        uncompressed[t][i].d[j] * uncompressed[t][i].t[j];


  double old_total_err = FLT_MAX;
  for(int it=0;it<num_it;it++)
  {
    memset(codebook_new, 0, sizeof(codebook_new));
    memset(cnt_new, 0, sizeof(cnt_new));
    double total_err = 0.0;
    // precompute partial smd of current codebook
    float codebook_mean[64] = {0};
    for(int i=0;i<64;i++) for(int j=0;j<8;j++)
      codebook_mean[i] += codebook_old[i].d[j] + 
        codebook_old[i].d[j] * codebook_old[i].t[j];

    // sort codebook by mean of codewords (use bubble-sort, we don't expect
    // things to change much after the first iteration):
    for(int c=0;c<64-1;c++)
    { // sort ascending order, codebook_mean[0] == 0
      int swapped = 0;
      for(int d=0;d<64-c-1;d++) if(codebook_mean[d] > codebook_mean[d+1])
      { // swap d and d+1
        float m = codebook_mean[d];
        codebook_mean[d] = codebook_mean[d+1];
        codebook_mean[d+1] = m;
        codeword_t w = codebook_old[d];
        codebook_old[d] = codebook_old[d+1];
        codebook_old[d+1] = w;
        swapped = 1;
      }
      if(!swapped) break;
    }

    for(int t=0;t<VOL_MOTION_SAMPLES;t++)
    { // quantise 2x2x2 blocks [i] for all time slices [t]
      for(int i=0;i<64;i++)
      {
        int ref = 0; // default to empty block
        float best_err = FLT_MAX;
        { // binary search on codebook_mean[]s to get closest to mean of this voxel block mean[t][i]
          int ref_m = 0, ref_M = 63;
          while(ref_M-ref_m > 1)
          {
            int m = (ref_M+ref_m)/2;
            const float err = codebook_mean[m] - mean[t][i];
            if(err < 0.0f) ref_m = m;
            else ref_M = m;
          }
          ref = ref_m;
          best_err = (codebook_mean[ref] - mean[t][i]) * (codebook_mean[ref] - mean[t][i]);
        }
        int go[2] = {1, 1};
        int pos[2] = {ref, ref};
        while(go[0] || go[1])
        {
          for(int dir=0;dir<2;dir++) if(go[dir])
          {
            pos[dir] += dir ? 1 : -1;
            if(pos[dir] >= 64 || pos[dir] < 0) go[dir] = 0;
            else
            {
              // calculate smd distance to codeword m
              const float d_mean = (mean[t][i] - codebook_mean[pos[dir]]) * (mean[t][i] - codebook_mean[pos[dir]]);
              if(d_mean > 16 * best_err) go[dir]= 0;
              else
              { // compute new best err, write ref if smaller than current min
                float err = 0.0f;
#ifdef _OPENMP
#pragma omp simd reduction(+:err)
#endif
                for(int j=0;j<8;j++)
                {
                  const float dens = uncompressed[t][i].d[j];
                  const float temp = uncompressed[t][i].t[j];
                  const float old_dens = codebook_old[pos[dir]].d[j];
                  const float old_temp = codebook_old[pos[dir]].t[j];
                  const float ds = dens - old_dens;
                  const float dt = dens * temp - old_dens * old_temp;
                  err += ds*ds + dt*dt;
                }
                if(err < best_err)
                {
                  ref = pos[dir];
                  best_err = err;
                }
              }
            }
          }
        }
        total_err += best_err;
        // accumulate our data to new codebook (weighted average)
        _vol_payload_compress_entry(p, i, t, ref);
        if(ref)
        { // do not average the zero block, we want it clean
#ifdef _OPENMP
#pragma omp simd
#endif
          for(int j=0;j<8;j++)
          {
            codebook_new[ref].d[j] += uncompressed[t][i].d[j];
            codebook_new[ref].t[j] += uncompressed[t][i].t[j];
          }
          cnt_new[ref]++;
        }
      }
    }
#if 0
    fprintf(stderr, "total error it %d : %g, counts:\n", it, total_err);
    for(int k=0;k<64;k++) fprintf(stderr, "%d ", cnt_new[k]);
    fprintf(stderr, "\n");
#endif
    for(int i=0;i<64;i++) if(cnt_new[i])
#ifdef _OPENMP
#pragma omp simd
#endif
      for(int k=0;k<8;k++)
      {
        codebook_old[i].d[k] = codebook_new[i].d[k]/cnt_new[i];
        codebook_old[i].t[k] = codebook_new[i].t[k]/cnt_new[i];
      }
    // happens for very boring blocks (mostly empty)
    if(total_err < 1e-6) break;
    if((old_total_err-total_err)/total_err < 0.04) break; // almost converged.
    old_total_err = total_err;
  }
  // write codebook_old to payload struct
  for(int i=0;i<512;i++)
  {
    const int i64 = ((i&0x6)>>1)|((i&0x30)>>2)|((i&0x180)>>3); // address of desired 4x4x4 block
    const int i8 = ((i&0x40)>>4)|((i&8)>>2)|(i&1);             // offset in 2x2x2 block
    p->d[i] = float_to_half(CLAMP(codebook_old[i64].d[i8], 0, 65504.0));
    p->t[i] = float_to_half(CLAMP(codebook_old[i64].t[i8], 0, 65504.0));
  }
  return isstatic;
}

// uncompress desired entry
static inline void vol_payload_uncompress(
    const vol_payload_compressed_t *p,
    int i, // spatial  index, 0<=i<512
    int t, // temporal index, 0<=t<VOL_MOTION_SAMPLES
    float *result)
{
  int i64 = ((i&0x6)>>1)|((i&0x30)>>2)|((i&0x180)>>3); // address of desired 4x4x4 block
  int off = i & 0x49; // offset in 2x2x2 block (in 512-base)

  int r;
#ifdef LOW_BITS
  int k = i64 * 6;         // bit address of 6-bit reference
  int kb = k / 24;         // block index (one block has 3 bytes)
  int ki  = (k - 24*kb)/6; // bit address within one block / 6 bit per number
  switch(ki)
  {
    case 0: // [0]:6
      r = p->ref[t][kb][0] & 0x3f;
      break;
    case 1: // [0]:2 and [1]:4
      r = (p->ref[t][kb][0] >> 6) | ((p->ref[t][kb][1] & 0xf) << 2);
      break;
    case 2: // [1]:4 and [2]:2
      r = ((p->ref[t][kb][1] & 0xf0) >> 4) | ((p->ref[t][kb][2] & 3) << 4);
      break;
    default: // case 3: [2]:6
      r = p->ref[t][kb][2] >> 2;
      break;
  }
#else
  r = p->ref[t][i64];
#endif

  int r512 = off + (((r&0x3)<<1)|((r&0xc)<<2)|((r&0x30)<<3)); // convert to 8x8x8 base address + offset
  result[0] = half_to_float(p->d[r512]);
  result[1] = half_to_float(p->t[r512]);
}
