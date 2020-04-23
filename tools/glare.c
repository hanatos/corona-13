#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdint.h>
// preempt any inclusion of colorspace related stuff:
#define CLAMP(a, m, M) fminf(M, fmaxf(m, a))
#include "../include/colorout_srgb.h"
#include "../include/colorin_ciergb.h"
#include "../include/spectrum.h"

#define M_PI           3.14159265358979323846  /* pi */
#define dot(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define cross(v1, v2, res) \
  (res)[0] = (v1)[1]*(v2)[2] - (v2)[1]*(v1)[2];\
  (res)[1] = (v1)[2]*(v2)[0] - (v2)[2]*(v1)[0];\
  (res)[2] = (v1)[0]*(v2)[1] - (v2)[0]*(v1)[1]

// units in cm
#define SENSOR_WD 1.0
#define SENSOR_HT 0.5
#define BOX_DIST 10.0

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

float sample_light(const float r, const float xp, float *x)
{
  x[0] = xp * 10.0;
  x[1] = 0.0f;
  x[2] = 10.0*BOX_DIST;
  return 1.0f;
}

#if 1 // box around stereo rig
float sample_edge(float r, float *x, float *e, int side)
{
  // position on the edge
  x[0] = (side ? 1.0 : -1.0) * (SENSOR_WD/2.0f);
  x[1] = (.5f-r) * SENSOR_HT;
  x[2] = BOX_DIST;

  // edge tangent. small perturbation with noise!
  e[0] = 0.0f;// + sinf(100.0f*x[1])*0.05f + sinf(1997*x[1])*0.03f;
  e[1] = 1.0f;// + cosf(137.0f*x[1])*0.05f;
  e[2] = 0.0f;
  const float ilen = 1.0f/sqrtf(fmaxf(0.0f, dot(e, e)));
  for(int k=0;k<3;k++) e[k] *= ilen;
  return 1.0f/SENSOR_HT;
}
#else // round lens housing
float sample_edge(float r, float *x, float *e)
{
  // position on the edge, sample hemicircle:
  const float phi = M_PI * r;
  x[0] =  sinf(phi) * SENSOR_WD/2.0;
  x[1] = -cosf(phi) * SENSOR_WD/2.0;
  x[2] = BOX_DIST;

  // edge tangent. TODO: small perturbation with noise!
  e[0] = cosf(phi);
  e[1] = sinf(phi);
  e[2] = 0.0f;
  return 1.0f/(M_PI*SENSOR_WD);
}
#endif

complex float sample_cone(const float r, const float *l, const float *x, const float *e, const float lambda, const int side, float *out)
{
  float in[3] = {x[0] - l[0], x[1] - l[1], x[2] - l[2]};
  const float inl2 = dot(in, in);
  const float inl = sqrtf(inl2);
  for(int k=0;k<3;k++) in[k] /= inl;
  const float cos_beta = dot(in, e);
  float normal[3], t[3];
  normal[0] = 0.0f;
  normal[1] = 0.0f;
  normal[2] = 1.0f; // vary this with noise, too? or not?
  cross(normal, e, t);

  // need to make sure t points towards the sensor:
  if(t[0] > 0.0f) for(int i=0;i<3;i++) t[i] = - t[i];

  // only sample lower hemisphere which points towards sensor.
  // delta is angle between out projected to n x t and n
  const float delta = (side ? 1.0 : 2.0) * M_PI/2.0f + r * M_PI/2.0f;
  // gamma is angle between in projected to n x t and n
  float inp[3] = {in[0], in[1], in[2]};
  float dot_in_e = dot(in, e);
  for(int k=0;k<3;k++) inp[k] -= dot_in_e * e[k];
  float ilen_inp = 1.0f/sqrtf(fmaxf(0.0f, dot(inp, inp)));
  for(int k=0;k<3;k++) inp[k] *= ilen_inp;
  // const float cos_gamma = fminf(1.0f, fmaxf(-1.0f, -dot(inp, normal)));
  // const float gamma = acosf(cos_gamma);

  const float sin_beta = sqrtf(fmaxf(0.0f, 1.0f-cos_beta*cos_beta));
  // geometric term connecting l and x:
  const float G = sin_beta/inl2;

  assert(t[0] <= 0.0f);
  for(int k=0;k<3;k++)
    out[k] = e[k] * cos_beta + sin_beta*(cosf(delta)*normal[k] + sinf(delta)*t[k]);

  // const float n = 0.5f; // 90 degree opening wedge, alpha = 3/2pi = pi(2-n) or n = 2 - alpha/pi
  // const float n = 1.75f; // pi/4
  // const float n = 1.95f;

#if 1
  // TODO: importance sample this stupid projected half angle 
  float outp[3] = {out[0], out[1], out[2]};
  const float dot_out_e = dot(out, e);
  for(int k=0;k<3;k++) outp[k] -= dot_out_e * e[k];
  const float ilen_outp = 1.0f/fmaxf(0.0f, sqrtf(dot(outp, outp)));
  for(int k=0;k<3;k++) outp[k] *= ilen_outp;
  const float cos_delta_gamma = CLAMP(-dot(inp, outp), -1.0, 1.0);
  const float delta_gamma = acosf(cos_delta_gamma);
#endif
  // use real keller formula here, also track phase:
  const complex float D = - cexpf(I * M_PI/4.0f)/(2.0f * sqrtf(1.0f/lambda) * sin_beta) *
    // + if boundary u = 0, - if boundary du/dn = 0
    (1.0f/cosf(delta_gamma/2.0f) + 1.0f/sinf(delta_gamma/2.0f));
  return 2.0f/M_PI * G * D;
#if 0
  assert(delta-gamma >= 0.0f);
  assert(delta-gamma <= M_PI);
  const float num = cosf(M_PI/4.0f) * sinf(M_PI/n);
  // const float den = n * M_PI * sin_beta;// * (cosf(M_PI/n) - cosf((delta-gamma)/n));
  const float den = n * M_PI * sin_beta * (cosf(M_PI/n) - cosf(delta_gamma/n));
  // 1/pdf(omega) * geometry term from light * term from keller's paper
  return 2.0f/M_PI * G * (num*num)/(den*den);
#endif
}

int main(int argc, char *argv[])
{
  const size_t width = 200, height = 100;
  float *out = (float*)malloc(sizeof(float)*3*width*height);
  complex float *buf = (complex float *)malloc(sizeof(complex float)*3*width*height);

  const size_t spp = 10000000;
  const int num_frames = 1;
  const int start_frame = 0;
  for(int f=start_frame;f<start_frame + num_frames;f++)
  {
    memset(buf, 0x0, sizeof(float)*3*width*height);
    float xp = 1.0 + f;
    size_t cnt = 0;
#pragma omp parallel for default(shared)
    for(size_t k=0;k<spp*width*height;k++)
    {
#pragma omp atomic
      cnt ++;
      if(!(cnt & 0x1fffff) || (cnt == spp*width*height-1)) fprintf(stderr, "\rrendering frame %d/%d %5.2f%%", f-start_frame+1, num_frames, 100.0*cnt/(spp*width*height));
      for(int side=0;side<2;side++)
      {
        float r0, r1, r2, r3;
        r0 = tea2(k, 0, &r1);
        r2 = tea2(k, 1, &r3);
        double lambda = spectrum_sample_lambda(r0, 0) * 1e-7; // cm
        float l[3], x[3], e[3], y[3];
        const float pdf_l = sample_light(r1, xp, l);
        const float pdf_x = sample_edge(r2, x, e, side);
        complex double throughput = width*height;//1e20; // light brightness
        throughput /= pdf_l;
        throughput /= pdf_x;
        // assert(throughput > 0.0f);

        float out[3];
        throughput *= sample_cone(r3, l, x, e, lambda, side, out);
        // assert(throughput > 0.0f);
        // trace ray (x,out) -> intersect sensor plane
        if(out[2] >= 0.0f) continue; // pointing away
        const double t = - x[2] / out[2];
        for(int i=0;i<3;i++) y[i] = x[i] + t * out[i];

        // apply phase shift:
        // 1e-7: no diff, 7e-6: maybe one mode, 1e-5: gone again
        throughput *= cexp(I*2.0*M_PI/lambda * t);
        // throughput *= cexp(I*2.0*M_PI/lambda *8e-6* t);
        // fprintf(stderr, "period %g\n", 2.0*M_PI/lambda * y[1]*1e-5 * 3e-2);//t);
        // throughput = t;

        // normalise by sample count
        throughput /= spp;
        const float pi = (y[0] + SENSOR_WD*.5f)/SENSOR_WD * width;
        const float pj = (y[1] + SENSOR_HT*.5f)/SENSOR_HT * height;
        const int pii = (int)pi;
        const int pji = (int)pj;

        // fprintf(stderr, "throughput %g\n", throughput);

        if(pii >= 0 && pii < width && pji >= 0 && pji < height)
        {
          float acc[6], rgb[3];
          spectrum_p_to_rgb(lambda*1e7, creal(throughput), rgb);
          acc[0] = rgb[0]; acc[2] = rgb[1], acc[4] = rgb[2];
          spectrum_p_to_rgb(lambda*1e7, cimag(throughput), rgb);
          acc[1] = rgb[0]; acc[3] = rgb[1], acc[5] = rgb[2];
          filter_accum(pi, pj, acc, (float *)buf, 6, 0, 6, width, height);
        }
      }
    }

    // convert to rgb:
    for(int k=0;k<3*width*height;k++)
      out[k] = cabsf(buf[k]);

    // write output image:
    char filename[1024];
    snprintf(filename, 1024, "glare_%04d.pfm", f);
    FILE *fp = fopen(filename, "wb");
    fprintf(fp, "PF\n%zu %zu\n-1.0\n", width, height);
    fwrite(out, sizeof(float)*3, width*height, fp);
    fclose(fp);
  }
  fprintf(stderr, "\n");
  free(buf);
  free(out);
  exit(0);
}
