#pragma once
#include "corona_common.h"
#include "points.h"

typedef enum vol_interpolation_t
{
  s_vol_constant    = 0<<0,
  s_vol_smooth      = 1<<0,
  s_vol_motion_blur = 1<<1,

  // keep nazicompiler happy:
  s_vol_constant_blur = s_vol_constant | s_vol_motion_blur,
  s_vol_smooth_blur = s_vol_smooth | s_vol_motion_blur,
}
vol_interpolation_t;

static inline void _vol_interpolate_sample(float *out)
{
  // Generation of Stratified Samples for B-Spline Pixel Filtering
  // Michael Stark Peter Shirley Michael Ashikhmin
  // unstratified recurrence / iterated box filter:
  const int n = 3; // cubic b-spline
  for(int k=0;k<4;k++)
  {
    out[k] =  - n/2.0;
    for(int i=0;i<=n;i++)
      out[k] += points_rand(rt.points, common_get_threadid());
  }
}

static inline void vol_interpolate_smooth(float *i, float *j, float *k, float *t, float sx, float st)
{
  float g[4];
  _vol_interpolate_sample(g);
  *i += g[0]*sx;
  *j += g[1]*sx;
  *k += g[2]*sx;
  *t  = fmodf(*t + g[3]*st, 1.0f);
}
