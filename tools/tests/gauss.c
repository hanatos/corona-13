#define DREGGN\
gcc -lm -D_GNU_SOURCE -Wall -I../../include gauss.c -o gauss; exit
#include "sampler_common.h"
#include "matrix2.h"
#include <string.h>

int main(int argc, char* argv[])
{
  // covariance matrix
  float S[4] = {.05f, .03f, .03f, .05f};
  float Si[4];
  mat2_invert(S, Si);

  float B[4] = {0.0f};
  B[0] = sqrtf(Si[0]);
  B[1] = Si[1]/B[0];
  B[2] = 0.0f;
  B[3] = sqrtf(Si[3] - B[1]*B[1]);
  fprintf(stderr, "B = %g %g %g %g\n", B[0], B[1], B[2], B[3]);
  float T[4], BT[4];
  mat2_transpose(B, BT);
  mat2_mul(BT, B, T);
  fprintf(stderr, "BT B = %g %g %g %g\n", T[0], T[1], T[2], T[3]);
  fprintf(stderr, "     = %g %g %g %g\n", Si[0], Si[1], Si[2], Si[3]);

  mat2_invert(B, T);

  const int wd = 512, ht = 512;
  float *pixel = (float *)malloc(sizeof(float)*wd*ht);
  memset(pixel, 0, sizeof(float)*wd*ht);

  const int N = 100*wd*ht;
  for(int k=0;k<N;k++)
  {
    float g[2], g2[2];
    float r1 = drand48(), r2 = drand48();
    sample_gaussian(r1, r2, g, g+1);
    mat2_mulv(T, g, g2);
    int i = (.5f + g2[0]) * wd, j = (.5f + g2[1]) * ht;
    if(i >= 0 && i < wd && j >= 0 && j < ht)
      pixel[i+wd*j] += wd*ht/(float)N;
  }
  // write
  FILE *f = fopen("sample.pfm", "wb");
  fprintf(f, "PF\n%d %d\n-1.0\n", wd, ht);
  for(int k=0;k<wd*ht;k++) for(int i=0;i<3;i++) fwrite(pixel+k, sizeof(float), 1, f);
  fclose(f);
  memset(pixel, 0, sizeof(float)*wd*ht);

  for(int j=0;j<ht;j++)
  {
    for(int i=0;i<wd;i++)
    {
      float vec[2] = {(i-wd/2.)/(float)wd, (j-ht/2.)/(float)ht};
      float tmp[2];
      mat2_mulv(Si, vec, tmp);
      float d = vec[0]*tmp[0] + vec[1]*tmp[1];
      pixel[i+wd*j] = 1.0f/(2.0f*M_PI * sqrtf(mat2_det(S))) * expf( - 0.5f * d);
    }
  }
  // write
  f = fopen("eval.pfm", "wb");
  fprintf(f, "PF\n%d %d\n-1.0\n", wd, ht);
  for(int k=0;k<wd*ht;k++) for(int i=0;i<3;i++) fwrite(pixel+k, sizeof(float), 1, f);
  fclose(f);
  free(pixel);
  exit(0);
}
