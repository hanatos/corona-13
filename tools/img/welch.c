#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#define MAX_P 1.f


// computes area left of t for t-distribution
// with v degrees of freedom and t negative
// from stack overflow:
static inline double tcdf64(const double t, const int v)
{
  // alg 3 appl. statistics 1968 vol 17 p 189
  const double x = v / (v + t * t);
  double c = 1.0f, s = 1.0f;
  int ioe = v % 2, k = 2 + ioe;
  if (v < 1) return 0.0f;
  if (v >= 4)
  {
    while(k <= v - 2)
    {
      c *= x - x / k;
      s += c;
      k += 2;
    }
  }
  c = t / sqrtf(v);
  if(ioe != 1)
    return 0.5f + 0.5f * sqrtf(x) * c * s;
  return 0.5f + ((v == 1 ? 0 : x * c * s) + atanf(c))/M_PI;
}

// from stack overflow:
static inline float tcdf(float t, int v)
{
  // alg 3 appl. statistics 1968 vol 17 p 189
  const float b = v / (v + t * t);
  float c = 1.0f, s = 1.0f;
  int ioe = v % 2, k = 2 + ioe;
  if (v < 1) return 0.0f;
  if (v >= 4)
  {
    while(k <= v - 2)
    {
      c *= b - b / k;
      s += c;
      k += 2;
    }
  }
  c = t / sqrtf(v);
  if(ioe != 1)
    return 0.5f + 0.5f * sqrtf(b) * c * s;
  return 0.5f + ((v == 1 ? 0 : b * c * s) + atanf(c))/M_PI;
}

static inline void viridis_quintic(
    float x,        // input float between 0 and 1
    float *col)     // output colour ramp
{
  x = fminf(1.0f, x);
  const float x2 = x*x;
  const float x3 = x2*x;
  const float x4 = x2*x2;
  const float x5 = x3*x2;
  col[0] = +0.280268003 -0.143510503 * x +2.225793877*x2 -14.815088879 * x3 +25.212752309 * x4 -11.772589584 * x5;
  col[1] = -0.002117546 +1.617109353 * x -1.909305070*x2 +2.701152864 *x3 -1.685288385 * x4 +0.178738871 *x5;
  col[2] = +0.300805501 +2.614650302 * x -12.019139090*x2 +28.933559110 *x3 -33.491294770* x4 +13.762053843 * x5;
}


float erf_approx3(float t)
{
  const float x = 1.0f / (1.0f + 0.4747f * t);
  const float a1 = 0.3480242;
  const float    a2 = -0.0958798;
  const float    a3 = 0.7478556;
  printf("t = %f, ", t);
  printf("x = %f, ", x);
  printf("e^(-x*x) = %f, ", exp(-x*x));
  printf("a1 + a2 + a3 = %f ", a1 + a2 + a3);
  const float val= 1.0 - (a1 * x + a2 * x * x + a3 * x * x * x) * exp(- t * t);
  printf("erf(%f)=%f, ", t, val);
  return val;
}


float phi(float x)
{
    return 0.5 * (1.f + erf_approx3(x / sqrt(2.f)));
}

float gauss_cdf_test(float x)
{
  if (x >= 0.f)
  {
    return 2.f * (1.f - phi(x));
  }
  return 2.f * (1.f - phi(-x));
}

static inline float gauss_cdf(float t)
{
  float tt = t > 0 ? t : -t;
  tt = tt / sqrt(2.0f);
  
  const float x = 1.f / (1.f + 0.47047 * tt);
  const float erf = 1.0 - x * (0.3480242 
                   + x*(-0.0958798
                        + 0.7478556 *x))
                 * exp(-tt*tt);
  return 1.f - erf;
}



void print_help()
{
    printf("Usage: ./welch firstimage.pfm secondimage.pfm <arguments>\n\n");

    printf("-p, --pscale <f>         scale p values by this\n");
    printf("-g                       replace t-dist with gaussian for integration\n");
    printf("-gt <i>                  set threshold to only replace t-dist with gaussian when nu >= this threshold\n");
    printf("--sqrt                   visualize sqrt of p-value instead of p-value\n");
    printf("-o <s>                   prefix of outfiles (will write out welch.pfm, pvalues.txt and nuvalues.txt\n");
    printf("--block-size <i>         length and with in pixels of pixel block used for one welch sample\n");
    printf("--help                   print this\n");
    printf("-3                       visualize all three color channels' p-values\n");
}

/*
 * === Usage:
 *
 * Pass first and second image as arguments
 *
 * If the image is called image_fb00.pfm there
 * should also be a file called image_welch (no extension).
 *
 * Also works with double frame buffer image_fbd00.pfm.
 *
 * The welch file doesn't need a header, the width
 * and height are extracted from the .pfm image.
 *
 *
 * === Neat resources:
 *
 * Sample Variance:
 * http://www.statisticshowto.com/sample-variance/
 *
 * Welch test intro and example:
 * https://www.youtube.com/watch?v=2-ecXltt2vI
 * https://www.youtube.com/watch?v=gzrmHpA54Sc
 *
 * p value table
 * https://math.stackexchange.com/questions/808474/p-value-for-lower-upper-tailed-t-test
 */
int main(int argc, char *argv[])
{

  float p_scale = 1.f;
  char * outname_prefix = "";
  int replace_with_gauss = 0;
  int gauss_nu_threshold = 0;
  int visualization_mode = 0; //0 = regular or scaled color scale, 1 = square root
  int block_size = 32;
  int vis_three = 0;

  if (argc < 2)
  {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++)
  {
    if ((strcmp(argv[i], "--pscale") == 0) && argc > i+1) p_scale = atof(argv[++i]);
    if ((strcmp(argv[i], "-p") == 0) && argc > i+1) p_scale = atof(argv[++i]);
    if ((strcmp(argv[i], "-g") == 0) ) replace_with_gauss = 1;
    if ((strcmp(argv[i], "-gt")==0) && argc > i+1) gauss_nu_threshold = atol(argv[++i]);
    if ((strcmp(argv[i], "--sqrt")==0)) visualization_mode = 1;
    if ((strcmp(argv[i], "-o") == 0) && argc > i+1) outname_prefix = argv[++i];
    if ((strcmp(argv[i], "--block-size")==0) && argc > i+1) block_size = atol(argv[++i]);
    if ((strcmp(argv[i], "--help")==0)) { print_help(); return 0; }
    if ((strcmp(argv[i], "-3")==0)) vis_three = 1;
  }
  
 // const float p_scale = (argc > 3) ? atof(argv[3]) : 1.f;
  printf("pscale %f\n",p_scale);
  printf("blocksize %d\n", block_size);
  printf("replace gauss? %d, threshold %d \n", replace_with_gauss, gauss_nu_threshold);
  printf("Write to prefix %s\n", outname_prefix);
  printf("Vis mode: %d\n", visualization_mode);
//  const char * outname_prefix = (argc > 4) ? argv[4] : "";
  

//  const int replace_with_gauss = (argc > 5) ? 1 : 0;
//  const int gauss_nu_threshold = 20; // only replace with gaussian if there are at least this many degrees of freedom.


  // read img1 and img2 (not needed anymore, except for width, height and debugging info)
  uint64_t width, height, width2, height2;
  FILE *fin = fopen(argv[1], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[1]); exit(1); }
  fscanf(fin, "PF\n%ld %ld\n%*[^\n]", &width, &height);
  fgetc(fin); // \n
  float *pixels1 = malloc(width*height*3*sizeof(float));
  fread(pixels1, sizeof(float), width*height*3, fin);
  fclose(fin);

  fin = fopen(argv[2], "rb");
  if(!fin) { fprintf(stderr, "could not open second %s!\n", argv[2]); exit(2); }
  fscanf(fin, "PF\n%ld %ld\n%*[^\n]", &width2, &height2);
  if(width2 != width || height2 != height) { fprintf(stderr, "resolution mismatch! height1: %lu, height2: %lu, width1: %lu, width2: %lu\n",
      height, height2, width, width2); exit(3); }
  fgetc(fin); // \n
  float *pixels2 = malloc(width*height*3*sizeof(float));
  fread(pixels2, sizeof(float), width*height*3, fin);
  fclose(fin);

  float *output = malloc(width*height*3*sizeof(float));
//  const int block_size = 32;

  const uint64_t w_wd = width / block_size;
  const uint64_t w_ht = height / block_size;

  // read sum of squares 1
  char filename[1024];
  strncpy(filename, argv[1], sizeof(filename));
  char *c = filename;
  while(strncmp(c, "_fb", 3) && c < filename+strlen(filename)) c++;
  if(c >= filename+strlen(filename)) { fprintf(stderr, "could not find welch buffers!\n"); exit(4);}
  sprintf(c, "_welch2");
  fin = fopen(filename, "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", filename); exit(5); }
  double *welch1_2 = malloc((w_wd*w_ht*3 + 2)*sizeof(double)); // last two values are sample count and overlays (samples per pixel; multiple of what counts as sample here)
  fread(welch1_2,sizeof(double), w_wd*w_ht*3 + 2, fin); // will gracefully crash if resolution doesn't match.
  fclose(fin);

  // read sum of squares 2
  strncpy(filename, argv[2], sizeof(filename));
  c = filename;
  while(strncmp(c, "_fb", 3) && c < filename+strlen(filename)) c++;
  if(c >= filename+strlen(filename)) { fprintf(stderr, "could not find welch buffers!\n"); exit(6);}
  sprintf(c, "_welch2");
  fin = fopen(filename, "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", filename); exit(7); }
  double *welch2_2 = malloc((w_wd*w_ht*3+2)*sizeof(double));
  fread(welch2_2, sizeof(double), w_wd*w_ht*3 + 2, fin); // will gracefully crash if resolution doesn't match.
  fclose(fin);


  // read sample sum 1
  strncpy(filename, argv[1], sizeof(filename));
  c = filename;
  while(strncmp(c, "_fb", 3) && c < filename+strlen(filename)) c++;
  if(c >= filename+strlen(filename)) { fprintf(stderr, "could not find welch buffers!\n"); exit(4);}
  sprintf(c, "_welch1");
  fin = fopen(filename, "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", filename); exit(5); }
  double *welch1_1 = malloc((w_wd*w_ht*3 + 2)*sizeof(double));
  fread(welch1_1,sizeof(double), w_wd*w_ht*3 + 2, fin); // will gracefully crash if resolution doesn't match.
  fclose(fin);

  // read sample sum 2
  strncpy(filename, argv[2], sizeof(filename));
  c = filename;
  while(strncmp(c, "_fb", 3) && c < filename+strlen(filename)) c++;
  if(c >= filename+strlen(filename)) { fprintf(stderr, "could not find welch buffers!\n"); exit(6);}
  sprintf(c, "_welch1");
  fin = fopen(filename, "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", filename); exit(7); }
  double *welch2_1 = malloc((w_wd*w_ht*3+2)*sizeof(double));
  fread(welch2_1, sizeof(double), w_wd*w_ht*3 + 2, fin); // will gracefully crash if resolution doesn't match.
  fclose(fin);


  double *pvals = malloc(w_wd*w_ht*sizeof(double)*3);
  int pvalcnt = 0; // only fill p-values where they were actually computed

  int *nuvals = malloc(w_wd*w_ht*sizeof(int)*3);

  const double welchsamples1 = welch1_2[w_wd*w_ht*3];
  const double welchsamples2 = welch2_2[w_wd*w_ht*3];
//  const double overlays1 = welch1_2[w_wd*w_ht*3+1]; // = welchsamples1 * 3 + {0,1,2}
//  const double overlays2 = welch2_2[w_wd*w_ht*3+1];

  const double n1 = welchsamples1;
  const double n2 = welchsamples2;

  printf("spp1 = %f, spp2 = %f\n", welch1_2[w_wd*w_ht*3+1], welch2_2[w_wd*w_ht*3+1]);
  printf("N1 = %f, N2 = %f\n", n1, n2);

  const double n1inv = n1 > 0.0 ? 1.0/n1 : 0.0;
  const double n2inv = n2 > 0.0 ? 1.0/n2 : 0.0;

  for(int j=0;j<w_ht;j++)
  {
    for(int i=0;i<w_wd;i++)
    {
      double sumX1[3] = {0.0};
      double sumX2[3] = {0.0};
      for(int k=0;k<3;k++)
      {
        sumX1[k] = welch1_1[3*(j*w_wd+i)+k];
        sumX2[k] = welch2_1[3*(j*w_wd+i)+k];
      }

      double s1_2[3], s2_2[3], tmp[3], t[3];
      int nu[3];
      float cdf = 1.0;
      float p_values[3];

      for(int k=0;k<3;k++)
      {
        // unbiased sample variance
        s1_2[k] = (1.0/(n1-1.0)) * (welch1_2[3*(j*w_wd+i)+k] - (n1inv*sumX1[k])*sumX1[k]);
        s2_2[k] = (1.0/(n2-1.0)) * (welch2_2[3*(j*w_wd+i)+k] - (n2inv*sumX2[k])*sumX2[k]);

        assert(s1_2[k] >= 0.0);
        assert(s2_2[k] >= 0.0);

        tmp[k] = s1_2[k] * n1inv + s2_2[k] * n2inv;
        if(tmp[k] == 0 || !(tmp[k] < FLT_MAX))
        {
          cdf = 1.0f;
          continue;
        }
        assert(tmp[k] >= 0.0);
        assert(tmp[k] == tmp[k]);

        // test lower bound of cdf interval:
        // t = (X_1 - X_2) / sqrt (s_1²/N_1 + s_2²/N_2)
        t[k] = -fabs(sumX1[k]*n1inv - sumX2[k]*n2inv) / sqrtf(tmp[k]);
        if(!(t[k] == t[k])
          || t[k] == 0 || !(t[k] < FLT_MAX))
        {
          cdf = 1.0f;
          continue;
        }

        // TODO: if these degrees of freedom are > 100 or so go full gaussian!
        nu[k] = roundf(tmp[k]*tmp[k] / (
              s1_2[k]*s1_2[k] * n1inv * n1inv / (n1-1.0f) +
              s2_2[k]*s2_2[k] * n2inv * n2inv / (n2-1.0f)
              ));
        if (i==0 && j == 0)
        {
          printf("tmp[%d] = %f\n",k, tmp[k]);
          printf("s1_2[%d] = %f\n", k, s1_2[k]);
          printf("s2_2[%d] = %f\n", k, s2_2[k]);
          printf("n1inv = %f\n", n1inv);
          printf("n2inv = %f\n", n2inv);
          printf("nu[%d] = %d\n", k, nu[k]);
        }
        //        printf("v[%d] = %d\n", k, nu[k]);
        // now test the null-hypothesis that the two means are equal, using
        // a two-tailed test on the t-distribution:
        const float p_value = (replace_with_gauss && nu[k] >= gauss_nu_threshold) ? gauss_cdf(t[k]) : 2.f * tcdf64(t[k], nu[k]);
        cdf = fminf(cdf, p_value); // this cdf is the p value. The smallest p value is at the "worst" color channel
        p_values[k] = p_value;

        assert(p_value>= 0.f && p_value <= 1.f);

        // if we do this we also get a p-value of 0 for every index where no p value was actually computed, e.g. if t[k] or tmp[k] = 0:
        // This has the advantage of knowing which index corresponds to each pixel, but will skew the p-value histogram.
        // pvals[3*(j*w_wd+i)+k] = p_value;

        // keep track of how many pvalues were actually computed. These values are useful for histogram, but we can't reconstruct
        // which value corresponds to which pixel.
        pvals[pvalcnt] = p_value;
        pvalcnt++;

        nuvals[3*(j*w_wd+i)+k] = nu[k];

      }
      if (vis_three)
      {
        for(int kk=0;kk<3;kk++)
        {
          float col[3];
          if (visualization_mode == 1)
            p_values[kk] = sqrt(p_values[kk]);
          const float confidence = fmaxf(0.f,1.f - (p_values[kk] > 0.f ? (p_scale*p_values[kk]) : 0.f));
          viridis_quintic(confidence/MAX_P, col);

          for(int jj=0;jj<block_size;jj++)
            for(int ii=(kk*block_size)/3;ii<((kk+1)*block_size)/3;ii++)
              for(int k=0;k<3;k++)
                output[3*((block_size*j+jj)*width + block_size*i+ii)+k] = col[k];
        }
      }
      else
      {
        float col[3];
        //the color method is yellow for large values and pink for values close to zero. The p value is good when it is close to one, so use 1-cdf as argument here
        if (visualization_mode == 1)
          cdf = sqrtf(cdf);
        const float confidence = fmaxf(0.f,1.f - (cdf > 0.f ? (p_scale*cdf) : 0.f));
        viridis_quintic(confidence/MAX_P, col);
        for(int jj=0;jj<block_size;jj++)
          for(int ii=0;ii<block_size;ii++)
            for(int k=0;k<3;k++)
              output[3*((block_size*j+jj)*width+block_size*i+ii)+k] = col[k];
      }
    }
  }

  // plot scale:
  for(int j=0;j<height;j++) for(int i=width-16;i<width-8;i++)
  {
    float f = 1.0f-j/(height-1.0f);
    for(int k=0;k<3;k++) output[3*(j*width+i)+k] = f;
  }
  for(int j=0;j<height;j++) for(int i=width-8;i<width;i++)
  {
    float f = 1.0f-j/(height-1.0f), col[3];
    viridis_quintic(f/MAX_P, col);
    for(int k=0;k<3;k++) output[3*(j*width+i)+k] = col[k];
  }


  char fnbuffer[1024];
  snprintf(fnbuffer, sizeof(fnbuffer), "%s_welch.pfm", outname_prefix);
  printf(fnbuffer);
  FILE *f = fopen(fnbuffer, "wb");
  if(f)
  {
    fprintf(f, "PF\n%ld %ld\n-1.0\n", width, height);
    fwrite(output, sizeof(float)*3, width*height, f);
    fclose(f);
    printf("finished, wrote results to %s\n", fnbuffer);
  }
  else
    printf("failed writing results.\n");

  snprintf(fnbuffer, sizeof(fnbuffer), "%s_pvalues.txt", outname_prefix);
  printf(fnbuffer);
  FILE *f1 = fopen(fnbuffer, "wb");
  if (f1)
  {
    // print number of p-values that were actually computed
    fprintf(f1, "%d ", pvalcnt);
    for(int k=0; k<pvalcnt; k++)
    {
      fprintf(f1, "%f ", pvals[k]);
    }
    fclose(f1);
  }

  snprintf(fnbuffer, sizeof(fnbuffer), "%s_nuvalues.txt", outname_prefix);
  printf(fnbuffer);
  FILE *f2 = fopen(fnbuffer, "wb");
  if (f2)
  {
    for (int k = 0; k < w_wd*w_ht*3; k++)
    {
      fprintf(f2, "%d ", nuvals[k]);
    }
    fclose(f2);
  }

#if 0 // debug print the cdf function
  for(int k=0;k<100;k++)
  {
    float t = -4.0 + 8.0*k/100.0f;
    fprintf(stderr, "%g %g\n", t, tcdf(t, 7));
  }
#endif
}

