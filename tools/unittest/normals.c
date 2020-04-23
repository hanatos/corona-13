// gcc -O0 -g -std=c11 -D_GNU_SOURCE normals.c -o normals -lm  -I ../../include
#include "prims_geo.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main(int argc, char *argv[])
{
  const int num = 1000000;
  float err_max = 0.0f;
  float err_sum = 0.0f;
  for(int k=0;k<num;k++)
  {
    float n[3] = {2.0f*drand48()-1.0f, 2.0f*drand48()-1.0f, 2.0f*drand48()-1.0f};
    // float n[3] = {0, 1, 1};
    float n2[3];
    normalise(n);
    const uint32_t enc = _geo_encode_normal(n);
    _geo_decode_normal(enc, n2);
    // float d[3] = {n[0] - n2[0], n[1] - n2[1], n[2] - n2[2]};
    // const float e = sqrtf(dotproduct(d,d));
    const float e = fabsf(acosf(CLAMP(dotproduct(n,n2), -1, 1)));
    err_max = fmaxf(e, err_max);
    err_sum += e;
  }
  fprintf(stderr, "normal vectors error: max = %g avg = %g\n", err_max, err_sum/num);


  // float / half for uv texture coords. > 10 we're not really interested in that.
  err_max = 0.0f;
  err_sum = 0.0f;
  for(int k=0;k<num;k++)
  {
    float f = 20.0f*drand48()-10.0f;
    // float f = (int)(4096.0f*drand48()-2048.0f); // these are precise :)
    const uint16_t enc = float_to_half(f);
    const float f2 = half_to_float(enc);
    const float e = fabsf(f - f2);
    err_max = fmaxf(e, err_max);
    err_sum += e;
  }
  fprintf(stderr, "half float error: max = %g avg = %g\n", err_max, err_sum/num);


  // uvw texture coords for hair strands, all in [0, 1):
  err_max = 0.0f;
  err_sum = 0.0f;
  for(int k=0;k<num;k++)
  {
    float n[3] = {drand48(), drand48(), drand48()};
    float n2[3];
    const uint32_t enc = _geo_encode_uvw(n[0], n[1], n[2]);
    _geo_decode_uvw(enc, n2, n2+1, n2+2);
    float d[3] = {n[0] - n2[0], n[1] - n2[1], n[2] - n2[2]};
    const float e = sqrtf(dotproduct(d,d));
    err_max = fmaxf(e, err_max);
    err_sum += e;
  }
  fprintf(stderr, "uvw error: max = %g avg = %g\n", err_max, err_sum/num);
  exit(0);
}
