#pragma once

static inline float lens_aperture_area(const float radius, const int blades)
{
  const float tri = .5f*radius * radius * sinf(2.0f*M_PI/(float)blades);
  return blades * tri;
}

static inline void lens_aperture_sample(float *x, float *y, float r1, float r2, const float radius, const int blades)
{
  const int tri = (int)(r1*blades);
  // rescale:
  r1 = r1*blades - tri;

  // sample triangle:
  float a = sqrtf(r1);
  float b = (1.0f-r2)*a;
  float c = r2*a;

  float p1[2], p2[2];

  common_sincosf(2.0f*M_PI/blades * (tri+1), p1, p1+1);
  common_sincosf(2.0f*M_PI/blades * tri, p2, p2+1);

  *x = radius * (b * p1[1] + c * p2[1]);
  *y = radius * (b * p1[0] + c * p2[0]);
}

static inline int lens_aperture_clip(const float x, const float y, const float radius, const int blades)
{ 
  // early out
  if(x*x + y*y > radius*radius) return 0;
  float xx = radius; 
  float yy = 0.0f;
  for(int b=1;b<blades+1;b++)
  {      
    float tmpx, tmpy;
    common_sincosf(2.0f*(float)M_PI/blades * b, &tmpy, &tmpx);
    tmpx *= radius;
    tmpy *= radius;
    const float normalx = xx + tmpx;
    const float normaly = yy + tmpy;
    float dot0 = (normalx)*(x-xx) + (normaly)*(y-yy);
    if(dot0 > 0.0f) return 0;
    xx = tmpx;
    yy = tmpy;
  }
  return 1;
}
