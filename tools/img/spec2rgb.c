#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
double drand48();
void srand48(long int);

#include "rgbe.h"
#include "spectrum.h"

// srgb dreggn
const int rgb_lambda_num = 89;
const float rgb_lambda_min = 380.0f;
const float rgb_lambda_step = 5.0f;

const float rgb_red[] = {
1.5000E-03, 3.8000E-03, 8.9000E-03, 1.8800E-02, 3.5000E-02, 5.3100E-02, 7.0200E-02, 7.6300E-02, 7.4500E-02, 5.6100E-02, 3.2300E-02, 4.4000E-03, 4.7800E-02, 9.7000E-02, 1.5860E-01, 2.2350E-01, 2.8480E-01, 3.3460E-01, 3.7760E-01, 4.1360E-01, 4.3170E-01, 4.4520E-01, 4.3500E-01, 4.1400E-01, 3.6730E-01, 2.8450E-01, 1.8550E-01, 4.3500E-02, 1.2700E-01, 3.1290E-01, 5.3620E-01, 7.7220E-01, 1.0059E+00, 1.2710E+00, 1.5574E+00, 1.8465E+00, 2.1511E+00, 2.4250E+00, 2.6574E+00, 2.9151E+00, 3.0779E+00, 3.1613E+00, 3.1673E+00, 3.1048E+00, 2.9462E+00, 2.7194E+00, 2.4526E+00, 2.1700E+00, 1.8358E+00, 1.5179E+00, 1.2428E+00, 1.0070E+00, 7.8270E-01, 5.9340E-01, 4.4420E-01, 3.2830E-01, 2.3940E-01, 1.7220E-01, 1.2210E-01, 8.5300E-02, 5.8600E-02, 4.0800E-02, 2.8400E-02, 1.9700E-02, 1.3500E-02, 9.2400E-03, 6.3800E-03, 4.4100E-03, 3.0700E-03, 2.1400E-03, 1.4900E-03, 1.0500E-03, 7.3900E-04, 5.2300E-04, 3.7200E-04, 2.6500E-04, 1.9000E-04, 1.3600E-04, 9.8400E-05, 7.1300E-05, 5.1800E-05, 3.7700E-05, 2.7600E-05, 2.0300E-05, 1.4900E-05, 1.1000E-05, 8.1800E-06, 6.0900E-06, 4.5500E-06};
const float rgb_green[] = {
-4.0000E-04, -1.0000E-03, -2.5000E-03, -5.9000E-03, -1.1900E-02, -2.0100E-02, -2.8900E-02, -3.3800E-02, -3.4900E-02, -2.7600E-02, -1.6900E-02, 2.4000E-03, 2.8300E-02, 6.3600E-02, 1.0820E-01, 1.6170E-01, 2.2010E-01, 2.7960E-01, 3.4280E-01, 4.0860E-01, 4.7160E-01, 5.4910E-01, 6.2600E-01, 7.0970E-01, 7.9350E-01, 8.7150E-01, 9.4770E-01, 9.9450E-01, 1.0203E+00, 1.0375E+00, 1.0517E+00, 1.0390E+00, 1.0029E+00, 9.6980E-01, 9.1620E-01, 8.5710E-01, 7.8230E-01, 6.9530E-01, 5.9660E-01, 5.0630E-01, 4.2030E-01, 3.3600E-01, 2.5910E-01, 1.9170E-01, 1.3670E-01, 9.3800E-02, 6.1100E-02, 3.7100E-02, 2.1500E-02, 1.1200E-02, 4.4000E-03, 7.8000E-05, -1.3680E-03, -1.9880E-03, -2.1680E-03, -2.0060E-03, -1.6420E-03, -1.2720E-03, -9.4700E-04, -6.8300E-04, -4.7800E-04, -3.3700E-04, -2.3500E-04, -1.6300E-04, -1.1100E-04, -7.4800E-05, -5.0800E-05, -3.4400E-05, -2.3400E-05, -1.5900E-05, -1.0700E-05, -7.2300E-06, -4.8700E-06, -3.2900E-06, -2.2200E-06, -1.5000E-06, -1.0200E-06, -6.8800E-07, -4.6500E-07, -3.1200E-07, -2.0800E-07, -1.3700E-07, -8.8000E-08, -5.5300E-08, -3.3600E-08, -1.9600E-08, -1.0900E-08, -5.7000E-09, -2.7700E-09};
const float rgb_blue[] = {
6.2000E-03 , 1.6100E-02 , 4.0000E-02 , 9.0600E-02 , 1.8020E-01 , 3.0880E-01 , 4.6700E-01 , 6.1520E-01 , 7.6380E-01 , 8.7780E-01 , 9.7550E-01 , 1.0019E+00 , 9.9960E-01 , 9.1390E-01 , 8.2970E-01 , 7.4170E-01 , 6.1340E-01 , 4.7200E-01 , 3.4950E-01 , 2.5640E-01 , 1.8190E-01 , 1.3070E-01 , 9.1000E-02 , 5.8000E-02 , 3.5700E-02 , 2.0000E-02 , 9.5000E-03 , 7.0000E-04 , 4.3000E-03 , 6.4000E-03 , 8.2000E-03 , 9.4000E-03 , 9.7000E-03 , 9.7000E-03 , 9.3000E-03 , 8.7000E-03 , 8.0000E-03 , 7.3000E-03 , 6.3000E-03 , 5.3700E-03 , 4.4500E-03 , 3.5700E-03 , 2.7700E-03 , 2.0800E-03 , 1.5000E-03 , 1.0300E-03 , 6.8000E-04 , 4.4200E-04 , 2.7200E-04 , 1.4100E-04 , 5.4900E-05 , 2.2000E-06 , 2.3700E-05 , 2.8600E-05 , 2.6100E-05 , 2.2500E-05 , 1.8200E-05 , 1.3900E-05 , 1.0300E-05 , 7.3800E-06 , 5.2200E-06 , 3.6700E-06 , 2.5600E-06 , 1.7600E-06 , 1.2000E-06 , 8.1700E-07 , 5.5500E-07 , 3.7500E-07 , 2.5400E-07 , 1.7100E-07 , 1.1600E-07 , 7.8500E-08 , 5.3100E-08 , 3.6000E-08 , 2.4400E-08 , 1.6500E-08 , 1.1200E-08 , 7.5300E-09 , 5.0700E-09 , 3.4000E-09 , 2.2700E-09 , 1.5000E-09 , 9.8600E-10 , 6.3900E-10 , 4.0700E-10 , 2.5300E-10 , 1.5200E-10 , 8.6400E-11 , 4.4200E-11};

static inline float get_spectrum(const float *s, float lambda)
{
  int l = (lambda - rgb_lambda_min)/rgb_lambda_step;
  if(l < 0 || l >= rgb_lambda_num) return 0.0f;
  return s[l];
}

int main(int argc, char *arg[])
{
#if 0
  // another of these endless tests.
  const int N = 100000, M = 31;
  float sum[3] = {0,0,0}, tmp[3], res[3];
#define func sinf(1.51*(lambda-400.0)/300.0)
// #define func 1.0
  srand48(666);
  for(int k=0;k<N;k++)
  {
    const float lambda = spectrum_sample_lambda(drand48(), NULL);
    spectrum_p_to_xyz(lambda, func, tmp);
    for(int k=0;k<3;k++) sum[k] += tmp[k]/(float)N;
  }
  spectrum_xyz_to_rgb(sum, res);
    printf("col %f %f %f\n", sum[0], sum[1], sum[2]);
    printf("col %f %f %f\n", res[0], res[1], res[2]);
  for(int k=0;k<=M;k++)
  {
    const float lambda = 400.0f + k/(float)M*(700.0-400.0);
    float p = spectrum_rgb_to_p(lambda, res);
    printf("p(%f) = %f (%f)\n", lambda, p, func);
  }
  exit(0);
#undef func
#endif

  char outname[512];
  char basename[512];
  if(argc < 5)
  {
    fprintf(stderr, "%s: usage inbase lambda_min lambda_step lambda_num [out.hdr]\n", arg[0]);
    exit(1);
  }
  if (argc == 5)
  {
    strncpy(outname, arg[1], sizeof(outname));
    strncpy(basename, arg[1], sizeof(basename));
    char* extension = outname + strlen(outname);
    strcpy(extension, ".hdr");
  }
  else
  {
    strncpy(outname, arg[5], sizeof(outname));
  }
  const float lambda_min = atof(arg[2]);
  const float lambda_step  = atof(arg[3]);
  const float lambda_num   = atof(arg[4]);
  float *pixels = NULL, *rgb = NULL;
  float col[3];
  int width, height;
  char filename[1024];
  for(int bin=0;bin<lambda_num;bin++)
  {
    const float lambda = lambda_min + lambda_step*bin;
    snprintf(filename, 1024, "%s_l%d.hdr", basename, (int)(lambda));
    // printf("loading %s\n", filename);
    FILE *f = fopen(filename, "rb");
    if(!f)
    {
      fprintf(stderr, "could not read hdri file: %s!\n", filename);
      return 1;
    }
    RGBE_ReadHeader(f, &width, &height, NULL);
    // assume these stay constant over the hdr files:
    if(width != 2*height) { fprintf(stderr, "ERROR: width has to be 2*height!\n"); return 1; }
    if(!pixels)
    {
      pixels = (float *)malloc(sizeof(float)*3*width*height);
      rgb = (float *)malloc(sizeof(float)*3*width*height);
      for(int i=0;i<3*width*height;i++) rgb[i] = 0.0f;
    }
    RGBE_ReadPixels_RLE(f, pixels, width, height);
    fclose(f);
    for(int j=0;j<height;j++) for(int i=0;i<width;i++)
    {
      int index = 3*(i + width*j);
      // spectrum_p_to_rgb(lambda, pixels[index], col);
      spectrum_p_to_xyz(lambda, pixels[index], col);
      // col[0] = get_spectrum(rgb_red, lambda);
      // col[1] = get_spectrum(rgb_green, lambda);
      // col[2] = get_spectrum(rgb_blue, lambda);
      for(int k=0;k<3;k++) rgb[index+k] += col[k]*pixels[index]/(lambda_num);
    }
  }
  FILE *f = fopen(outname, "wb");
  if(f)
  {
    for(int i=0;i<3*width*height;i+=3)
    {
      spectrum_xyz_to_rgb(rgb+i, col);
      for(int k=0;k<3;k++) rgb[i+k] = col[k];
    }
    RGBE_WriteHeader(f, width, height, NULL);
    RGBE_WritePixels(f, rgb, width*height);
    fclose(f);
  }
  free(rgb);
  free(pixels);
  exit(0);
}
