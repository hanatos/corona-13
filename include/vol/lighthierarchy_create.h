#pragma once

#include "vol/types.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

typedef struct vol_payload_flux_t
{
  // store normalised importance sampling weights for flux per voxel of internal octree.
  // that is mu_e * L_e * volume of voxel, normalised to total sum in this octree block.
  // data layout is p[8] is the parent node, p[0-7] are the bottom level probabilities:
  // every node has a binary tree layout: [ z | y0 y1 | x00 x01 x10 x11 ]
  uint16_t p[9][7];
  uint16_t q[9][7];
}
vol_payload_flux_t;

static inline float _lighthierarchy_store(
    uint16_t *const buf16,  // output: 7 quantised values
    const float *const v)   // input values, 8 voxels
{
  for(int k=0;k<8;k++) assert(v[k] == v[k]);
  float binary[7] = {0.0f};
  float total = v[0]+v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7], t = 0.;
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
  for(int k=0;k<7;k++)
    buf16[k] = CLAMP(binary[k]*0x10000, 0, 0xffff);
  return total;
}

// fill one light hierarchy node
static inline void vol_lighthierarchy_create_node(
    vol_payload_uncompressed_t *parent, // parent block
    vol_payload_uncompressed_t *child,  // uncompressed child block to aggregate into parent voxel
    int idx,                            // index of the given child voxel inside the parent block
    vol_payload_flux_t *const node,     // output node
    int leaf_level,                     // flag whether we're called on leaf level (and need the shader)
    vol_emission_shader_t shader,       // shader to compute L_e
    int force_static)
{
  float fluxt_avg[8] = {0.0f};
  float fluxt_max[8] = {0.0f};
  if(parent) for(int t=0;t<VOL_MOTION_SAMPLES;t++)
  {
    parent->d[idx][t] = 0;
    parent->t[idx][t] = FLT_MAX; // for minorants
  }
  // this would only parallelise up to 8 threads, in the best case. not a good idea:
// #ifdef _OPENMP
// #pragma omp parallel for collapse(3) schedule(dynamic) default(shared)
// #endif
  for(int tz=0;tz<2;tz++) for(int ty=0;ty<2;ty++) for(int tx=0;tx<2;tx++)
  {
    float fluxb_avg[8] = {0.0f};
    float fluxb_max[8] = {0.0f};
    for(int bz=0;bz<2;bz++) for(int by=0;by<2;by++) for(int bx=0;bx<2;bx++)
    {
      for(int fz=0;fz<2;fz++) for(int fy=0;fy<2;fy++) for(int fx=0;fx<2;fx++)
      {
        int cidx = (4*tz+2*bz+fz)*64 + (4*ty+2*by+fy)*8 + (4*tx+2*bx+fx);
        float curr_flux_avg = 0.0f;
        float curr_flux_max = 0.0f;
        const int mcnt = (force_static?1:VOL_MOTION_SAMPLES);
        for(int t=0;t<mcnt;t++)
        { // accumulate over time
          if(child->d[cidx][t] > 65504.0f)// || (child->d[cidx][t] > 0.0f && child->d[cidx][t] < 1e-5f))
            fprintf(stderr, "[lh create] d out of half range: %g\n", child->d[cidx][t]);
          if(child->t[cidx][t] > 65504.0f)// || (child->t[cidx][t] > 0.0f && child->t[cidx][t] < 1e-5f))
            fprintf(stderr, "[lh create] t out of half range: %g\n", child->t[cidx][t]);
          const float density     = half_to_float(float_to_half(CLAMP(child->d[cidx][t], 0, 65504)));
          const float temperature = half_to_float(float_to_half(CLAMP(child->t[cidx][t], 0, 65504)));
          float flux_avg = temperature;
          float flux_max = density;
          if(leaf_level)
          {
            flux_avg = flux_max = 0.0f;
            // TODO: maybe have to run shaders on whole block and in SSE?
            if(density > 0.0 && temperature > 0.0) for(int l=0;l<16/MF_COUNT;l++)
            { // average over wavelengths
              float lf[MF_COUNT];
              for(int l=0;l<MF_COUNT;l++)
                lf[l] = fmodf(drand48() + l/(float)MF_COUNT, 1.0f);
              const mf_t lambda = spectrum_sample_lambda(mf_loadu(lf), 0);
              flux_max += density * mf_hsum(shader(density, temperature, lambda))/16.0f;
            }
            flux_avg = flux_max;
            // if(parent) parent->d[idx][t] = MAX(parent->d[idx][t], flux_avg);
            // if(parent) parent->t[idx][t] += flux_avg / (512.0f * mcnt);
            if(parent) parent->d[idx][t] = MAX(parent->d[idx][t], density); // majorants for woodcock
            if(parent) parent->t[idx][t] = MIN(parent->t[idx][t], density); // minorant for residual tracking
          }
          else if(parent)
          {
            // parent->d[idx][t] += density     / (512.0f * mcnt);
            // parent->t[idx][t] += temperature / (512.0f * mcnt);
            parent->d[idx][t] = MAX(parent->d[idx][t], density); // majorants for woodcock/ratio tracking
            parent->t[idx][t] = MIN(parent->t[idx][t], density); // minorants for residual tracking
          }
          curr_flux_avg += flux_avg;
          if (leaf_level)
            curr_flux_max = MAX(curr_flux_max, flux_max);
          else
            curr_flux_max += flux_max;
        }
        fluxb_avg[4*bz+2*by+bx] += curr_flux_avg;
        if (leaf_level)
          fluxb_max[4*bz+2*by+bx] = MAX(fluxb_max[4*bz+2*by+bx], curr_flux_max);
        else
          fluxb_max[4*bz+2*by+bx] += curr_flux_max;
      }
    }
    fluxt_avg[4*tz+2*ty+tx] = _lighthierarchy_store(node->p[4*tz+2*ty+tx], fluxb_avg);
    fluxt_max[4*tz+2*ty+tx] = _lighthierarchy_store(node->q[4*tz+2*ty+tx], fluxb_max);
  }
  _lighthierarchy_store(node->p[8], fluxt_avg);
  _lighthierarchy_store(node->q[8], fluxt_max);
}

