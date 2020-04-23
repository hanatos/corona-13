#pragma once
#include "prims.h"
#include <stdio.h>

// wraps a texture for easier passing around:
typedef struct texture_t
{
  float *tex;
  int channels;
  int width, height;

  float *mip_buf;
  float **mip;
  int mip_levels;
}
texture_t;

static inline void texture_lookup_mip(
    const texture_t *t,          // texture struct
    float u0, float v0,          // uv coordinates for full res
    int level,                   // mipmap level, 0 is full res image
    float *res)
{
  level = CLAMP(level, 0, t->mip_levels-1);
  const float mul = 1.0f/(1<<level);
  u0 = (u0 - floorf(u0))*t->width;
  v0 = (v0 - floorf(v0))*t->height;
  int u = u0 * mul + .5f, v = v0 * mul + .5f;
  const int w = t->width >> level, h = t->height >> level;
  float *buf = t->tex;
  if(level) buf = t->mip[level];
  const int c = t->channels;
  if(u >=0 && u < w && v >= 0 && v < h)
    for(int i=0;i<c;i++)
      res[i] = buf[c*(w*(h-v-1) + u)+i];
}

static inline void texture_lookup_mip_bilinear(
    const texture_t *t,          // texture struct
    float u0, float v0,          // uv coordinates for full res
    int level,                   // mipmap level, 0 is full res image
    float *res)
{
  level = CLAMP(level, 0, t->mip_levels-1);
  const float mul = 1.0f/(1<<level);
  u0 = (u0 - floorf(u0))*t->width;
  v0 = (v0 - floorf(v0))*t->height;
  float u = u0 * mul, v = v0 * mul;
  const int w = t->width >> level, h = t->height >> level;
  float *buf = t->tex;
  if(level) buf = t->mip[level];
  const int c = t->channels;
  const int iu = (int)u;
  const int iv = (int)v;
  const float fracu = u - iu, fracv = v - iv;
  if(iu >=0 && iu < w-1 && iv >= 0 && iv < h-1)
    for(int i=0;i<c;i++)
      res[i] =
        buf[c*(w*(h-iv-0-1) + iu + 0)+i] * (1.0f-fracu)*(1.0f-fracv) +
        buf[c*(w*(h-iv-0-1) + iu + 1)+i] * (     fracu)*(1.0f-fracv) +
        buf[c*(w*(h-iv-1-1) + iu + 0)+i] * (1.0f-fracu)*(     fracv) +
        buf[c*(w*(h-iv-1-1) + iu + 1)+i] * (     fracu)*(     fracv);
}

static inline void texture_lookup_bilinear(
    const texture_t *t,          // texture struct
    float u, float v,            // uv coordinates
    float *res)
{
  texture_lookup_mip_bilinear(t, u, v, 0, res);
}

static inline void texture_lookup_trilinear(
    const texture_t *t,
    float u, float v,
    float level,
    float *res)
{
  assert(t->mip);
  assert(t->mip_levels);

  int l0 = (int)level, l1 = (int)(level+1);
  int lf = level - l0;

  float r0[t->channels], r1[t->channels];
  texture_lookup_mip_bilinear(t, u, v, l0, r0);
  texture_lookup_mip_bilinear(t, u, v, l1, r1);
  for(int k=0;k<t->channels;k++) res[k] = (1.0f-lf) * r0[k] + lf * r1[k];
}

static inline void texture_lookup_ewa(
    const texture_t *t,    // texture struct
    float u, float v,      // uv coordinates
    const float d0[2],     // main axis of elliptical lookup
    const float d1[2],     // second axis of elliptical lookup
    float *res)
{
  float tmp[t->channels];
  for(int i=0;i<t->channels;i++) res[i] = 0.0f;
  float p[18] = { -1.0,  0.0,
                   0.0,  1.0,
                   1.0,  0.0,
                   0.0, -1.0,
                   0.0,  0.0,
                  -0.48,  0.52,
                   0.52,  0.48,
                   0.48, -0.52,
                  -0.52, -0.48 };
  float w[9] = {1/16., 1/16., 1/16., 1/16., 4/16., 2/16., 2/16., 2/16., 2/16.};
  const float max_sidelen = MAX(fabsf(d0[0]*t->width), fabsf(d0[1]*t->height));
  // const float max_sidelen = MIN(fabsf(d1[0]*t->width), fabsf(d1[1]*t->height));
  if((int)max_sidelen == 0)
  { // magnification
    texture_lookup_mip(t, u, v, 0, res); // nearest
    // float p0[2], p1[2];
    // p0[0] = 1.0/t->width;
    // p0[1] = 0.0;
    // p1[0] = 0.0;
    // p1[1] = 1.0/t->height;
    // for(int i=0;i<t->channels;i++) res[i] = 0.0f;
    // for(int k=0;k<18;k+=2)
    // {
    //   texture_lookup_bilinear(t, u+p0[0]*p[k]+p1[0]*p[k+1], v+p0[1]*p[k]+p1[1]*p[k+1], tmp);
    //   for(int i=0;i<t->channels;i++) res[i] += w[k/2] * tmp[i];
    // }
  }
  else
  { // minification
    // select mip map level based on shorter axis, to avoid overblur:
    const uint32_t sidelen_px = (int)MAX(fabsf(d1[0]*t->width), fabsf(d1[1]*t->height));
    const float l = 32 - __builtin_clz(sidelen_px);
    for(int i=0;i<t->channels;i++) res[i] = 0.0f;
    for(int k=0;k<18;k+=2)
    {
      texture_lookup_trilinear(t, u+d0[0]*p[k]+d1[0]*p[k+1], v+d0[1]*p[k]+d1[1]*p[k+1], l, tmp);
      for(int i=0;i<t->channels;i++) res[i] += w[k/2] * tmp[i];
    }
  }
}

static inline void texture_create_mip(
    texture_t *t)
{
  int num_levels = 0;
  for(int m = MIN(t->width, t->height); m; m>>=1) num_levels++;
  t->mip_levels = num_levels;
  t->mip = (float **)malloc(sizeof(float *) * num_levels);
  t->mip[0] = t->tex;
  size_t mem = 0;
  for(int l=1, w=t->width/2, h=t->height/2;l<num_levels;l++,w >>= 1,h >>= 1)
    mem += w*h * sizeof(float)*t->channels;
  t->mip_buf = (float *)malloc(mem);

  size_t off = 0;
  for(int l=1, w=t->width/2, h=t->height/2;l<num_levels;l++,w >>= 1,h >>= 1)
  {
    t->mip[l] = t->mip_buf + off;
    off += w*h * t->channels;

    // fill pixel data. should really use `Non-Power-of-Two Mipmap Creation'
    // http://http.download.nvidia.com/developer/Papers/2005/NP2_Mipmapping/NP2_Mipmap_Creation.pdf
    for(int j=0;j<h;j++) for(int i=0;i<w;i++)
      for(int c=0;c<t->channels;c++)
        t->mip[l][t->channels*(w*j+i)+c] = .25 * (
            t->mip[l-1][t->channels*((t->width>>(l-1))*(2*j+0)+2*i+0)+c] +
            t->mip[l-1][t->channels*((t->width>>(l-1))*(2*j+0)+2*i+1)+c] +
            t->mip[l-1][t->channels*((t->width>>(l-1))*(2*j+1)+2*i+0)+c] +
            t->mip[l-1][t->channels*((t->width>>(l-1))*(2*j+1)+2*i+1)+c]);
#if 0
    char fn[256];
    snprintf(fn, 256, "mip_%d.pfm", l);
    FILE *o = fopen(fn, "wb");
    if(t->channels == 3) fprintf(o, "PF\n%d %d\n-1.0\n", w, h);
    else fprintf(o, "Pf\n%d %d\n-1.0\n", w, h);
    fwrite(t->mip[l], sizeof(float)*t->channels, w*h, o);
    fclose(o);
#endif
  }
}

static inline int texture_load(
    texture_t *t,                // your struct, possibly on stack
    const char *filename,        // texture file name
    const int mip)               // create mip map levels or not
{
  memset(t, 0, sizeof(texture_t));
  FILE *f = fopen(filename, "rb");
  if(!f) return 1;
  if(fgetc(f) != 'P') goto error;
  const char type = fgetc(f);
  if(type == 'f') t->channels = 1; // grey level
  else if(type == 'F') t->channels = 3; // colour
  else if(type == 'x') t->channels = 1; // will read actual count later
  else goto error;

  if(type == 'x')
  {
    if(fscanf(f, "\n%d %d\n%d", &t->width, &t->height, &t->channels) != 3) goto error;
    // fprintf(stderr, "[texture] loading %dx%dx%d texture\n", t->width, t->height, t->channels);
  }
  else
  {
    if(fscanf(f, "\n%d %d\n%*[^\n]", &t->width, &t->height) != 2) goto error;
  }
  fgetc(f); // munch new line
  t->tex = (float *)malloc(t->channels*sizeof(float)*t->width*t->height);
  if(fread(t->tex, t->channels*sizeof(float), t->width*t->height, f) != t->width*t->height) goto error;
  fclose(f);

  if(mip) texture_create_mip(t);
  return 0;
error:
  free(t->tex);
  fclose(f);
  return 2;
}

static inline int texture_write(
    const texture_t *t,
    const char *filename)
{
  FILE *f = fopen(filename, "wb");
  if(!f) return 1;
  fprintf(f, "Px\n%d %d\n%d\n", t->width, t->height, t->channels);
  fwrite(t->tex, t->channels*sizeof(float), t->width*t->height, f);
  fclose(f);
  return 0;
}

static inline void texture_alloc(
    texture_t *t,
    const int width,
    const int height,
    const int channels)
{
  memset(t, 0, sizeof(texture_t));
  t->width = width;
  t->height = height;
  t->channels = channels;
  t->tex = (float *)malloc(sizeof(float)*width*height*channels);
}

static inline void texture_cleanup(
    texture_t *t)
{
  free(t->tex);
  free(t->mip);
  free(t->mip_buf);
}
