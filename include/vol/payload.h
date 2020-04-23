#pragma once

#include "vol/shaders.h"
#include "spectrum.h"
#include "half.h"
#include "corona_common.h" // for CLAMP
#include <math.h>
#include <assert.h>
#include <float.h>

// voxel volumes: all payload specific stuff should go into here.
// this is decoupled from the tree structure so we can experiment with it more easily


#define VOL_MOTION_SAMPLES 64
// #define VOL_MOTION_SAMPLES 4
// #define VOL_MOTION_SAMPLES 8


typedef struct vol_payload_uncompressed_t
{
  float t[512][VOL_MOTION_SAMPLES]; // temperature
  float d[512][VOL_MOTION_SAMPLES]; // density
}
vol_payload_uncompressed_t;

#include "vol/payload_compress.h"

typedef enum vol_payload_type_t
{
  s_vol_empty = 0,
  s_vol_static = 1,
  s_vol_full = 2,
}
vol_payload_type_t;

// client function to fill voxel i in payload data, using world space coordinate p[3]:
typedef vol_payload_type_t (* vol_payload_fill_t)(void *data, vol_payload_uncompressed_t *payload, const float *aabb, const int force_static);

static inline void vol_payload_get(
    const vol_payload_t *payload,  // data block to fill
    const int idx,                 // address of voxel
    const float time,              // ray time [0,1]
    const int isstatic,            // flag whether payload is static
    float *result)                 // write data here
{
  if(isstatic)
  {
    result[0] = half_to_float(payload->d[idx]);
    result[1] = half_to_float(payload->t[idx]);
    return;
  }
  const int bin = CLAMP(time * VOL_MOTION_SAMPLES, 0, VOL_MOTION_SAMPLES-1);
  vol_payload_uncompress(payload, idx, bin, result);
}

