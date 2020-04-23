/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "denoise_wavelets.h"

/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */


typedef float elem_type ;

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }


/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


elem_type kth_smallest(elem_type a[], int n, int k)
{
    register int i,j,l,m ;
    register elem_type x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}


#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


int main(int argc, char **arg)
{
  char outname[512], varname[512];
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s in.pfm var.pfm [out.ppm]\n", arg[0]);
    return 1;
  }
  if (argc < 4)
  {
    if (strlen(arg[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, arg[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, "-denoised.ppm");
    strncpy(varname, arg[1], sizeof(varname));
    extension = varname + strlen(varname) - 4;
    strcpy(extension, "-var.pfm");
  }
  else
  {
    strncpy(varname, arg[2], sizeof(varname));
    strncpy(outname, arg[3], sizeof(outname));
  }

  FILE *in = fopen(arg[1], "rb");
  FILE *out = fopen(outname, "wb");
  FILE *var = fopen(varname, "rb");
  if(!in || !out || !var)
  {
    fprintf(stderr, "could not open in and out files!\n");
    return 2;
  }

  int w, h;
  fscanf(in, "PF\n%d %d\n%*[^\n]\n", &w, &h);
  fscanf(var, "PF\n%d %d\n%*[^\n]\n", &w, &h);
  float *image  = (float *)malloc(sizeof(float)*3*w*h);
  float *varimg = (float *)malloc(sizeof(float)*3*w*h);
  unsigned char *img = (unsigned char *)malloc(sizeof(unsigned char)*3*w*h);
  fread(image,  sizeof(float)*3*w*h, 1, in);
  fread(varimg, sizeof(float)*3*w*h, 1, var);

  // transform to yuv
  for(int k=0;k<w*h;k++)
  {
    float rgb[3] = {image[3*k], image[3*k+1], image[3*k+2]};
    image[3*k+0] =  0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2];
    image[3*k+1] = -0.147*rgb[0] - 0.289*rgb[1] + 0.437*rgb[2];
    image[3*k+2] =  0.615*rgb[0] - 0.515*rgb[1] - 0.100*rgb[2];
  }

  // wavelet transform input, invariance
  const int numl = 6;
  const float prior[][3] = {{0,0,0},
                           // {6.253728, 2.263844, 5.404001},
                           // {4.920455, 1.583165, 3.160959},
                           // {2.902243, 0.517697, 0.799946},
                           // {3.794065, 0.122645, 0.176320}};
{10.751644, 2.360202, 5.999412},
{7.794510,  1.739266, 3.589897},
{3.775210,  0.560354, 0.888908}};

  float **tmp = (float **)malloc(sizeof(float *)*numl);
  for(int k=1;k<numl;k++)
  {
    const int wd = (int)(1 + (w>>(k-1))), ht = (int)(1 + (h>>(k-1)));
    tmp[k] = (float *)malloc(sizeof(float)*wd*ht);
  }

  // for(int k=0;k<3*w*h;k++) image[k] = varimg[k];
  for(int level=1;level<numl;level++) dt_iop_equalizer_wtf(image,  tmp, level, w, h, 1);
  // for(int level=1;level<numl;level++) dt_iop_equalizer_wtf(varimg, tmp, level, w, h, 0);

  // print out statistics
  int cnt[numl];
  float band_hist[numl][3];
  for(int i=0;i<numl;i++) cnt[i] = 0;
  for(int l=1;l<numl;l++)
  {
    band_hist[l][0] = band_hist[l][1] = band_hist[l][2] = 0.0f;
    cnt[l]++;
    for(int ch=0;ch<3;ch++)
    {
      const int step = 1<<l;
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) band_hist[l][ch] += image[3*w*j + 3*i + ch]*image[3*w*j + 3*i + ch];
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      band_hist[l][ch] += image[3*w*j + 3*i + ch]*image[3*w*j + 3*i + ch];
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) band_hist[l][ch] += image[3*w*j + 3*i + ch]*image[3*w*j + 3*i + ch];
    }
  }
  float band_max = 0.0f;
  for(int i=0;i<numl;i++) for(int ch=0;ch<3;ch++)
  {
    if(cnt[i]) band_hist[i][ch] /= cnt[i];
    else band_hist[i][ch] = 0.0;
    band_max = fmaxf(band_max, band_hist[i][ch]);
  }

  for(int ch=0;ch<3;ch++)
  {
    printf("channel %d\n", ch);
    for(int l=1;l<numl;l++)
    {
      printf("%f ", band_hist[l][ch]);
    }
    printf("\n");
  }

  // estimate sigma by Median(|Yij|)/0.6745 in HH1 (finest details)
  float noise_sigma[3];
  float noise_gain[3] = {3.5, 4.0, 4.0}; // set your 1/(noise tolerance) here
  // float noise_gain[3] = {3.5, 10.0, 10.0}; // set your 1/(noise tolerance) here
  {
    const int l = 1;
    const int step = 1<<l;
    float *hh1 = (float *)malloc(sizeof(float)*w/2*h/2);
    for(int ch=0;ch<3;ch++)
    {
      int ind = 0;
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) hh1[ind++] = fabsf(image[3*w*j + 3*i + ch]);
      noise_sigma[ch] = median(hh1, w/2*h/2)/0.6745;
      noise_sigma[ch] *= noise_gain[ch];
      printf("noise sigma channel %d = %f\n", ch, noise_sigma[ch]);
    }
  }


  const float blend = 0.;
  // attenuate wavelet details based on transformed variance
  for(int l=1;l<numl;l++)
  {
    for(int ch=0;ch<3;ch++)
    {
      const int step = 1<<l;
#if 0
      // coefficients in range [0, 2], 1 being neutral.
#define coeff (ch>0?4.0:1.0)*sqrtf(varimg[3*w*j+3*i+0] + varimg[3*w*j+3*i+1] + varimg[3*w*j+3*i+2])/(step*step*step)
// #define coeff 3./(step)
      // for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] -= coeff;
      // for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      image[3*w*j + 3*i + ch] -= coeff;
      // for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] -= 2*coeff; //*coeff;
#undef coeff
#endif

      // estimate noise variance for this sub-band
      float sigma_y_2, sigma_x, T;
      int sig_cnt;

      sigma_y_2 = 0.0;
      sig_cnt = 0;
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) { sigma_y_2 += (image[3*w*j + 3*i + ch])*(image[3*w*j + 3*i + ch]); sig_cnt ++;}
      sigma_y_2 /= sig_cnt;
      sigma_x = sqrtf(fmaxf(sigma_y_2 - noise_sigma[ch]*noise_sigma[ch], 0.0f));
      T = noise_sigma[ch]*noise_sigma[ch]/sigma_x;
      printf("LH T = %f = %f / %f\n", T, noise_sigma[ch]*noise_sigma[ch], sigma_x);
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);

      sigma_y_2 = 0.0;
      sig_cnt = 0;
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)     { sigma_y_2 += (image[3*w*j + 3*i + ch])*(image[3*w*j + 3*i + ch]); sig_cnt ++;}
      sigma_y_2 /= sig_cnt;
      sigma_x = sqrtf(fmaxf(sigma_y_2 - noise_sigma[ch]*noise_sigma[ch], 0.0f));
      T = noise_sigma[ch]*noise_sigma[ch]/sigma_x;
      printf("HL T = %f\n", T);
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)     image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);

      sigma_y_2 = 0.0;
      sig_cnt = 0;
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) { sigma_y_2 += (image[3*w*j + 3*i + ch])*(image[3*w*j + 3*i + ch]); sig_cnt ++;}
      sigma_y_2 /= sig_cnt;
      sigma_x = sqrtf(fmaxf(sigma_y_2 - noise_sigma[ch]*noise_sigma[ch], 0.0f));
      T = noise_sigma[ch]*noise_sigma[ch]/sigma_x;
      printf("HH T = %f\n", T);
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);


#if 0
      float sigma = 0.0;
      int sig_cnt = 0;
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) { sigma += fabsf(varimg[3*w*j + 3*i + ch]); sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      { sigma += fabsf(varimg[3*w*j + 3*i + ch]); sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) { sigma += fabsf(varimg[3*w*j + 3*i + ch]); sig_cnt ++;}
      sigma /= sig_cnt;

      // estimate signal variance sigma_x
      float mu_x = 0.0;
      sig_cnt = 0;
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) { mu_x += image[3*w*j + 3*i + ch]; sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      { mu_x += image[3*w*j + 3*i + ch]; sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) { mu_x += image[3*w*j + 3*i + ch]; sig_cnt ++;}
      mu_x /= sig_cnt;
      float sigma_x = 0.0;
      sig_cnt = 0;
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) { sigma_x += (image[3*w*j + 3*i + ch]-mu_x)*(image[3*w*j + 3*i + ch]-mu_x); sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      { sigma_x += (image[3*w*j + 3*i + ch]-mu_x)*(image[3*w*j + 3*i + ch]-mu_x); sig_cnt ++;}
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) { sigma_x += (image[3*w*j + 3*i + ch]-mu_x)*(image[3*w*j + 3*i + ch]-mu_x); sig_cnt ++;}
      sigma_x /= sig_cnt;

      sigma_x = sqrtf(fmaxf(0.001, sigma_x - sigma));

      // bayes shrink threshold:
#define T (sigma / sigma_x)
// #define T (fabsf(varimg[3*w*j + 3*i + ch]) / sigma_x)
      printf("T = %f = %f/%f\n", T, sigma, sigma_x);
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*copysignf(fmaxf(0, fabsf(image[3*w*j + 3*i + ch]) - T), image[3*w*j + 3*i + ch]);
#undef T
#endif

#if 0
      // prior histogram matching:
      for(int j=0;j<h;j+=step)      for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*prior[l][ch]*image[3*w*j + 3*i + ch]/band_hist[l][ch];
      for(int j=step/2;j<h;j+=step) for(int i=0;i<w;i+=step)      image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*prior[l][ch]*image[3*w*j + 3*i + ch]/band_hist[l][ch];
      for(int j=step/2;j<h;j+=step) for(int i=step/2;i<w;i+=step) image[3*w*j + 3*i + ch] = blend*image[3*w*j + 3*i + ch] + (1.0 - blend)*prior[l][ch]*image[3*w*j + 3*i + ch]/band_hist[l][ch];
#endif
    }
  }

  // backtransform
  for(int level=numl-1;level>0;level--) dt_iop_equalizer_iwtf(image,  tmp, level, w, h);
  // for(int level=numl-1;level>0;level--) dt_iop_equalizer_iwtf(varimg, tmp, level, w, h);
 
  for(int k=0;k<w*h;k++)
  {
    float yuv[3] = {image[3*k], image[3*k+1], image[3*k+2]};
    image[3*k+0] = yuv[0]                + 1.140*yuv[2];
    image[3*k+1] = yuv[0] - 0.394*yuv[1] - 0.581*yuv[2];
    image[3*k+2] = yuv[0] + 2.028*yuv[1];
  }

  for(int k=1;k<numl;k++) free(tmp[k]);
  free(tmp);

  for(int k=0;k<3*w*h;k++) img[k]   = (unsigned char)(255.0f*(fminf(1.0, fmaxf(0.0, image[k]))));
  fprintf(out, "P6\n%d %d\n255\n", w, h);
  for(int j=h-1;j>=0;j--) fwrite(img + 3*w*j, 1, sizeof(unsigned char)*w*3, out);

  fclose(in);
  fclose(out);
  fclose(var);
  free(image);
  free(img);
  free(varimg);
  return 0;
}
