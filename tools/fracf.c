#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <stdint.h>
// #define MIN(a, b) ((a) > (b) ? (b) : (a))
// #define MAX(a, b) ((a) < (b) ? (b) : (a))
// #define CLAMP(a, m, M) MIN(MAX(m, a), M)
// #include "../include/colorin_ciergb.h"
// #include "../include/colorout_srgb.h"
#include "../include/spectrum.h"


float
tea2(unsigned int v0, unsigned int v1, float *out)
{
  const unsigned int key0 = 0xa341316c;
  const unsigned int key1 = 0xc8013ea4;
  const unsigned int key2 = 0xad90777d;
  const unsigned int key3 = 0x7e95761e;
  unsigned int sum = 0;
  unsigned int delta = 0x9e3779b9;

  // this can be lowered for performance,
  // 16 results in reasonably good points.
  const unsigned int rounds = 16;
  for(unsigned int i = 0; i < rounds; i++)
  {
    sum += delta;
    v0 += ((v1 << 4) + key0) ^ (v1 + sum) ^ ((v1 >> 5) + key1);
    v1 += ((v0 << 4) + key2) ^ (v0 + sum) ^ ((v0 >> 5) + key3);
  }
  // convert to float in [0, 1)
  *out = v1*(1.0f/4294967296.0f);
  return v0*(1.0f/4294967296.0f);
}

#if 0
// violates strict aliasing
static inline float tofloat(uint32_t i)
{
  return *(float *)&i;
}

static inline uint32_t touint(float f)
{
  return *(uint32_t *)&f;
}

static inline void common_atomic_add(float *f, float inc)
{
  uint32_t *i = (uint32_t *)f;
  float newval;
  uint32_t newi, oldi;
  do
  {
    newval = f[0];
    oldi = touint(newval);
    newval += inc;
    newi = touint(newval);
  }
  while(!__sync_bool_compare_and_swap(i, oldi, newi));
}
#endif

static inline float filter_bh_w(const float n)
{
  const float NN = 4.0f;
  if(n > NN-1.0f || n < 0.0f) return 0.0f;
  const float a0 = 0.35875;
  const float a1 = 0.48829;
  const float a2 = 0.14128;
  const float a3 = 0.01168;
  const float N_1 = 1.0f/(NN-1.0f);
  const float cos1 = cosf(2.0f*M_PI*n*N_1);
  const float cos2 = cosf(4.0f*M_PI*n*N_1);
  const float cos3 = cosf(6.0f*M_PI*n*N_1);
  return a0 - a1*cos1 + a2*cos2 - a3*cos3;
}

static inline void filter_accum(const float i, const float j, const float *const rad, float *pixel, const int bpp, const int offs, const int channels, const int width, const int height)
{
  int x = (int)i;
  int y = (int)j;
  float fx = i - x;
  float fy = j - y;

  const float dx = (int)(fx - 1.5f);
  const float dy = (int)(fy - 1.5f);

  const int x0 = (int)(i + dx);
  const int y0 = (int)(j + dy);
  const int u0 = x0 >= 0 ? 0 : - x0;
  const int v0 = y0 >= 0 ? 0 : - y0;
  const int u4 = x0 + 4 > width  ? width  - x0 : 4;
  const int v4 = y0 + 4 > height ? height - y0 : 4;
  for(int v=v0;v<v4;v++) for(int u=u0;u<u4;u++)
  {
    const float uu = u + fx - 1.f, vv = v + fy - 1.f;
    const float r = sqrtf(uu*uu + vv*vv);
    const float f = filter_bh_w(r+1.5f);
    for(int k=0;k<channels;k++)
#if 1//def FILTER_ATOMIC
      common_atomic_add(pixel+bpp*(x0+u+width*(y0+v)) + offs + k, rad[k]*f);
#else
    pixel[bpp*(x0+u+rt.width*(y0+v)) + offs + k] += rad[k]*f;
#endif
  }
}

int main(int argc, char *arg[])
{
  FILE *fp = fopen(arg[1], "rb");
  if(!fp)
  {
    fprintf(stderr, "[%s] couldn't open `%s'\n", arg[0], arg[1]);
    exit(1);
  }
  // init buf
  int width, height;
  int rd = fscanf(fp, "PF\n%d %d\n%*[^\n]", &width, &height);
  if(rd != 2)
  {
    fprintf(stderr, "[%s] `%s' doesn't look like a pfm image\n", arg[0], arg[1]);
    exit(2);
  }
  fgetc(fp);
  float *buf = (float *)aligned_alloc(128, sizeof(float)*width*height*3);
  rd = fread(buf, sizeof(float)*width*height*3, 1, fp);
  fclose(fp);

  complex float *f = (complex float *)aligned_alloc(128, sizeof(complex float)*width*height);
  complex float *F = (complex float *)aligned_alloc(128, sizeof(complex float)*width*height);

  const float D = 1.0f;//.05f; // distance, in micrometers
  const float lambda[3] = {0.670, 0.550, 0.400}; // wavelength in micrometers
  const complex float ipil[3] = { I*M_PI/(lambda[0] * D), I*M_PI/(lambda[1] * D), I*M_PI/(lambda[2] * D)};
  // this can apparently not go higher than 10x-20x D on a 230x230 image => aliasing.
  const float size = 10.0f; // in micrometers
  const float bandlimit = size*0.1f;
  const float bl2 = bandlimit*bandlimit;

  fprintf(stderr, "[%s] processing %d x %d image", arg[0], width, height);
  // TODO: sample lambda instead?
  // direct n^4 version:
  for(int c=0;c<3;c++)
  {
    // depends on distance to screen, far field would be pi/2 => fourier transform, fraunhofer diffraction
    // const float alpha = 0.001f; // pi/2 : fourier transform, multiple of pi: spatial domain
    // const float csca = 1.0f/sinf(alpha);
    // const float cota = 1.0f/tanf(alpha);

    fprintf(stderr, "\n[%s] processing color channel %d/3\n", arg[0], c+1);
    // init f from buf[c]
    for(int k=0;k<width*height;k++)
      f[k] = buf[3*k+c];

#if 1
    // fresnel diffraction, MC sampled
    size_t cnt = 0;
    // is about 100x less than 512^4
    const size_t num_samples = 1ul<<30;
    const float norm = 1.0;//width*height/num_samples;
    memset(F, 0, sizeof(complex float)*width*height);

#pragma omp parallel for schedule(static) default(shared)
    for(size_t k=0;k<num_samples;k++)
    {
      // const int tid = omp_get_thread_num();
      // uu, vv, xx, yy are all in micrometers
      float rand[2];
      rand[1] = tea2(k>>30, 2*(k&0xfffffff), rand);
      const float vf = rand[0]*height;
      const float uf = rand[1]*width;
      //onst int v = CLAMP((int)vf, 0, height-1);
      // const int u = CLAMP((int)uf, 0, width-1);
      const float vv = (rand[0] - .5f)*size;
      const float uu = (rand[1] - .5f)*size;
      rand[1] = tea2(k>>30, 2*(k&0xfffffff)+1, rand);
      // TODO: respect bounding box of aperture
      const float xm = 0.4f, xM = 0.6f, ym = .4f, yM = .6f;
      const int y = CLAMP((int)(rand[0]*height), 0, height-1);
      const int x = CLAMP((int)(rand[1]*width), 0, width-1);
      const float yy = (rand[0] - .5f)*size;
      const float xx = (rand[1] - .5f)*size;
#pragma omp atomic
      cnt++;
      if(!(cnt & 0xff) || (cnt == num_samples-1))
        fprintf(stderr, "\r[%s] processing %03.2f%%", arg[0], 100.0*cnt/((double)num_samples-1.0));

      // skip black pixels of input
      if(crealf(f[y*width+x]) <= 0.0) continue;

      // const complex float accum = norm * cexpf(-ipil[c] * (xx*xx + yy*yy) + 2.0f*ipil[c] * (uu*xx + vv*yy)) * f[y*width+x];
      const complex float chirp = cexpf(ipil[c] * (
                                         fminf(bl2, 2.0f*(uu*xx + vv*yy)) -
                                         fminf(bl2, (xx*xx + yy*yy))));
      const complex float accum = norm * chirp * f[y*width+x];
#if 1
      filter_accum(uf, vf, (float *)&accum, (float *)F, 2, 0, 2, width, height);
#else
#pragma omp critical
      {
        F[width*v + u] += accum;
      }
#endif
    }
    fprintf(stderr, "\n[%s] normalising", arg[0]);
    // normalisation:
#pragma omp parallel for schedule(static) default(shared)
    for(int k=0;k<width*height;k++)
    {
      const int v = k / width;
      const int u = k - v*width;
      // uu, vv, xx, yy are all in micrometers
      const float vv = (v/(height-1.0) - .5f)*size;
      const float uu = (u/(width-1.0) - .5f)*size;
      F[k] *= ipil[c]/M_PI * cexpf(-ipil[c] * fminf(bl2, (uu*uu+vv*vv)));
    }

    // init buf[c] from F
    for(int k=0;k<width*height;k++)
      buf[3*k+c] = cabsf(F[k]);
  }
#endif

#if 0
    const float norm = 1.0f/(width*height);
    int cnt = 0;

#pragma omp parallel for schedule(static) default(shared)
    for(int k=0;k<width*height;k++)
    {
      const int v = k / width;
      const int u = k - v*width;
      // uu, vv, xx, yy are all in micrometers
      const float vv = (v/(height-1.0) - .5f)*size;
#pragma omp atomic
      cnt++;
      if(!(cnt & 0xff) || (cnt == width*height-1))
        fprintf(stderr, "\r[%s] processing %03.2f%%", arg[0], 100.0f*cnt/(width*height-1.0f));
      const float uu = (u/(width-1.0) - .5f)*size;
      F[width*v + u] = 0.0f;
      for(int y=0;y<height;y++)
      {
        const float yy = (y/(height-1.0) - .5f)*size;
        for(int x=0;x<width;x++)
        {
          const float xx = (x/(width-1.0) - .5f)*size;
          // wiki, probably wrong:
          // F[width*v + u] += cexpf(-I*M_PI/lambda[c]*csca * (uu*xx + vv*yy) + I * cota*.5f/lambda[c] * (xx*xx + yy*yy)) * f[y*width + x];
          // fractional fourier transform from the paper:
          // F[width*v + u] += cexpf(I/lambda[c]*csca * (uu*xx + vv*yy)
                                 // -I/2.0 * cota/lambda[c] * (xx*xx + yy*yy)) * f[y*width + x];
          // fresnel diffraction:
          F[width*v + u] += norm * cexpf(-ipil[c] * (xx*xx + yy*yy) + 2.0f*ipil[c] * (uu*xx + vv*yy)) * f[y*width+x];
        }
      }
      // wiki:
      // F[width*v + u] *= (1.0f-I*cota) * cexpf(I*M_PI/lambda[c]*cota * (uu*uu + vv*vv));
      // paper:
      // F[width*v + u] *= I/lambda[c] * csca /(2.0f*M_PI) * cexpf(- I*0.5f*cota / lambda[c] * (uu*uu + vv*vv) - I*alpha);
      // fresnel diffraction:
      F[width*v + u] *= ipil[c]/M_PI * cexpf(-ipil[c] * (uu*uu+vv*vv));
    }

    // init buf[c] from F
    for(int k=0;k<width*height;k++)
      buf[3*k+c] = cabsf(F[k]);
  }
#endif
#if 0
    // factored version, looks different.
    // transform the image f(x,y) to fractional fourier domain F(u,v)
    const float norm = 1.0f;// /(width*height);
    // vertical, f -> F
    for(int v=0;v<height;v++)
    {
      fprintf(stderr, "\r[%s] processing row %d/%d", arg[0], v, height);
      for(int u=0;u<width;u++)
      {
        const float uu = (u/(width-1.0) - .5f)*size;
        F[width*v + u] = 0.0f;
        for(int x=0;x<width;x++)
        {
          const float xx = (x/(width-1.0) - .5f)*size;
          // F[width*v + u] += norm * cexpf(-ipil[c] * (xx*xx) + 2.0f*ipil[c] * (uu*xx)) * f[v*width+x];
          F[width*v + u] += norm * cexpf(ipil[c] * (xx*xx) - 2.0f*ipil[c] * (uu*xx)) * f[v*width+x];
          // F[width*v + u] += cexpf(-I*M_PI*csca * uu * xx + I * cota/2 * xx*xx) * f[v*width + x];
        }
        // F[width*v + u] *= csqrtf(1.0f-I*cota) * cexpf(I*M_PI*cota * uu*uu);
      }
    }

    // horizontal F -> f
    for(int u=0;u<width;u++)
    {
      fprintf(stderr, "\r[%s] processing col %d/%d", arg[0], u, width);
      for(int v=0;v<height;v++)
      {
        const float vv = (v/(height-1.0) - .5f)*size;
        f[width*v + u] = 0.0f;
        for(int y=0;y<height;y++)
        {
          const float yy = (y/(height-1.0) - .5f)*size;
          // f[width*v + u] += norm * cexpf(-ipil[c] * (yy*yy) + 2.0f*ipil[c] * (vv*yy)) * F[y*width+u];
          f[width*v + u] += norm * cexpf(ipil[c] * (yy*yy) - 2.0f*ipil[c] * (vv*yy)) * F[y*width+u];
          // f[width*v + u] += cexpf(-I*M_PI*csca * vv * yy + I * cota/2 * yy*yy) * F[y*width + u];
        }
        // f[width*v + u] *= csqrtf(1.0f-I*cota) * cexpf(I*M_PI*cota * vv*vv);
      }
    }

    // normalisation:
#pragma omp parallel for schedule(static) default(shared)
    for(int k=0;k<width*height;k++)
    {
      const int v = k / width;
      const int u = k - v*width;
      // uu, vv, xx, yy are all in micrometers
      const float vv = (v/(height-1.0) - .5f)*size;
      const float uu = (u/(width-1.0) - .5f)*size;
      f[k] *= ipil[c]/M_PI * cexpf(-ipil[c] * (uu*uu+vv*vv));
    }

    // init buf[c] from f
    for(int k=0;k<width*height;k++)
      buf[3*k+c] = cabsf(f[k]);
  }
#endif
  fprintf(stderr, "\n");

  // write buf
  fp = fopen("output.pfm", "wb");
  fprintf(fp, "PF\n%d %d\n-1.0\n", width, height);
  fwrite(buf, width*height, 3*sizeof(float), fp);
  fclose(fp);

  exit(0);
}
