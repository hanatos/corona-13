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
#include <assert.h>
void srand48(long);
double drand48();

#include "shader.h"
#include "rgbe.h"
#include "spectrum.h"
#include "sampler_common.h"

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

typedef struct envmap_t
{
  float *pixels;
  float **mip;
  double sum;
}
envmap_t;

typedef struct rgbe_t
{
  envmap_t *envmap;
  int levels;
  int width, w2n;
  int height, h2n;
  float aspectx, aspecty;
  float lambda_min, lambda_step;
  int lambda_num;
  float phase;
  int select;
}
rgbe_t;

int get_bin(rgbe_t *t, const float lambda)
{
  int l = (lambda - t->lambda_min)/t->lambda_step;
  if(l < 0) return 0;
  if(l >= t->lambda_num) return t->lambda_num - 1;
  return l;
}

float sky_rgbeis_sh(float *col)
{
  // return fmaxf(fmaxf(col[0], col[1]), col[2]);
  return col[0];// + col[1] + col[2];
}

void cleanup(void *data)
{
  rgbe_t *t = (rgbe_t *)data;
  for(int k=0;k<t->lambda_num;k++)
  {
    free(t->envmap[k].pixels);
    free(t->envmap[k].mip);
  }
  free(t->envmap);
  free(t);
}

float sky(rayhit_t *hit, float *dir, void* p)
{
  return 0.0f;
}

float add_sun(const rayhit_t *hit, const ray_t *ray, void *data)
{
  rgbe_t *t = (rgbe_t *)data;
  const int bin = get_bin(t, hit->lambda);
  // if(t->select && bin != t->select) return 0.0;
  float col[3];
  float x, y;
  if(fabsf(ray->dir[2]) > 0.999f)
  {
    y = 0;
    x = 0;
  }
  else
  {
    y = acosf(ray->dir[2])/M_PI * (t->height-1);
    x = fmodf(t->phase + M_PI + atan2f(ray->dir[0], ray->dir[1]), 2.0f*M_PI)/(2*M_PI) * (t->width-1);
  }
  if(!(x < t->width && x >= 0.0 && y < t->height && y >= 0.0)) x = y = 0.0;
  const int i = (int)x, j = (int)y;
#if 0
  col[0] = t->pixels[3*(t->width*j+i)];
  col[1] = t->pixels[3*(t->width*j+i)+1];
  col[2] = t->pixels[3*(t->width*j+i)+2];
#else
  const float u = x - i, v = y - j;
  const int jp = j == t->height - 1 ? j : j+1;
  const int ip = i == t->width  - 1 ? i : i+1;
  col[0] = (1-u)*((1-v)*t->envmap[bin].pixels[3*(t->width*j+i)] + v*t->envmap[bin].pixels[3*(t->width*jp+i)]) +
      u*((1-v)*t->envmap[bin].pixels[3*(t->width*j+ip)] + v*t->envmap[bin].pixels[3*(t->width*jp+ip)]);
  col[1] = (1-u)*((1-v)*t->envmap[bin].pixels[3*(t->width*j+i)+1] + v*t->envmap[bin].pixels[3*(t->width*jp+i)+1]) +
      u*((1-v)*t->envmap[bin].pixels[3*(t->width*j+ip)+1] + v*t->envmap[bin].pixels[3*(t->width*jp+ip)+1]);
  col[2] = (1-u)*((1-v)*t->envmap[bin].pixels[3*(t->width*j+i)+2] + v*t->envmap[bin].pixels[3*(t->width*jp+i)+2]) +
      u*((1-v)*t->envmap[bin].pixels[3*(t->width*j+ip)+2] + v*t->envmap[bin].pixels[3*(t->width*jp+ip)+2]);
#endif
  // compensate premultiplied sin
  // return (1.0f/sinf(M_PI*(.5f+y)/(float)t->height))*spectrum_rgb_to_p(hit->lambda, col);
  return (1.0f/sinf(M_PI*(.5f+y)/(float)t->height))*col[0];
}

float sample_sun(rayhit_t *hit, ray_t *ray, const float x1, const float x2, float *pdf, void *data)
{
  rgbe_t *t = (rgbe_t *)data;
  const int bin = get_bin(t, hit->lambda);
  // if(t->select && bin != t->select) return 0.0;
  double wx = x1, wy = x2;
#if 1
  if(wx < t->envmap[bin].mip[t->levels-1][0])
    wx *= .5/t->envmap[bin].mip[t->levels-1][0];
  else
    wx = 1.0 - (1.0-wx)*.5/(1.0 - t->envmap[bin].mip[t->levels-1][0]);
  wx *= 2.0;
  for(int m=t->levels-2;m>=0;m--)
  {
    int i = (int)wx, j = (int)wy, wd = t->w2n>>m;
    if(2*i >= wd) i = (wd-1)/2;
    if(4*j >= wd) j = (wd-1)/4;
    if(wy < j + (t->envmap[bin].mip[m][2*i + 2*j*wd]+t->envmap[bin].mip[m][2*i+1 + 2*j*wd]))
    {
      wy = j + (wy-j)*.5/(t->envmap[bin].mip[m][2*i + 2*j*wd]+t->envmap[bin].mip[m][2*i+1 + 2*j*wd]);
      j = 2*j;
    }
    else if (t->envmap[bin].mip[m][2*i + 2*j*wd]+t->envmap[bin].mip[m][2*i+1 + 2*j*wd] < 0.999)
    {
      wy = j + 1.0 - (1.0-(wy-j))*.5/(1.0 - (t->envmap[bin].mip[m][2*i + 2*j*wd]+t->envmap[bin].mip[m][2*i+1 + 2*j*wd]));
      //wy = j + 1.0 - (1.0-(wy-j))*.5/(t->mip[m][2*i + (2*j+1)*wd]+t->mip[m][2*i+1 + (2*j+1)*wd]);
      j = 2*j + 1;
    }
    const double f = t->envmap[bin].mip[m][2*i + j*wd]/(t->envmap[bin].mip[m][2*i + j*wd] + t->envmap[bin].mip[m][2*i+1 + j*wd]);
    if(wx <= i + f)
      wx = i + (wx-i)*.5/f;
    else if(wx > i + f) // do nothing for nan
      wx = i + 1.0 - (1.0-wx+i)*.5/(1.0-f);
    wx *= 2.0; wy *= 2.0f;
  }
#endif
  float col[3];
  const int i = (int)(wx*t->aspectx), j = (int)(wy*t->aspecty);
  col[0] = t->envmap[bin].pixels[3*(t->width*j+i)];
  col[1] = t->envmap[bin].pixels[3*(t->width*j+i)+1];
  col[2] = t->envmap[bin].pixels[3*(t->width*j+i)+2];
  const float theta = M_PI*(wy/(float)t->h2n);
  const float phi   = 2*M_PI*(wx/(float)t->w2n) - M_PI - t->phase; // - M_PI for imp sampling
  // printf("theta, phi %f %f\n", theta, phi);
  ray->dir[0] = sinf(phi)*sinf(theta);
  ray->dir[1] = cosf(phi)*sinf(theta);
  ray->dir[2] = cosf(theta);
  // if(pdf) *pdf = sky_rgbeis_sh(col)/(M_PI*M_PI*2.0*t->sum);
  // return M_PI*M_PI*2.0*t->envmap[bin].sum/*sinf(theta)*/*spectrum_rgb_to_p(hit->lambda, col)/sky_rgbeis_sh(col);
  if(pdf) *pdf = col[0]/(M_PI*M_PI*2.0*t->envmap[bin].sum);
  return M_PI*M_PI*2.0*t->envmap[bin].sum;
}

float sun_pdf(const rayhit_t *hit, const ray_t *ray, void *data)
{
  rgbe_t *t = (rgbe_t *)data;
  const int bin = get_bin(t, hit->lambda);
  // if(t->select && bin != t->select) return 0.0;
  float col[3];
  float x, y;
  if(fabsf(ray->dir[2]) > 0.999f)
  {
    y = 0;
    x = 0;
  }
  else
  {
    y = acosf(ray->dir[2])/M_PI * (t->height-1);
    x = fmodf(t->phase + M_PI + atan2f(ray->dir[0], ray->dir[1]), 2.0f*M_PI)/(2*M_PI) * (t->width-1);
  }
  if(!(x < t->width && x >= 0.0 && y < t->height && y >= 0.0)) x = y = 0.0;
  const int i = (int)x, j = (int)y;
  col[0] = t->envmap[bin].pixels[3*(t->width*j+i)];
  col[1] = t->envmap[bin].pixels[3*(t->width*j+i)+1];
  col[2] = t->envmap[bin].pixels[3*(t->width*j+i)+2];
  // return sky_rgbeis_sh(col)/(M_PI*M_PI*2.0*t->sum);
  return col[0]/(M_PI*M_PI*2.0*t->envmap[bin].sum); // constant grey for all channels.
}

int init(FILE* s, void **data)
{
  rgbe_t *t = (rgbe_t *)malloc(sizeof(rgbe_t));
  *data = t;
  char basename[1024], filename[1024];
  char line[1024];
  int dreggn = 0;
  dreggn = fscanf(s, "%[^\n]\n", line);
  float b = 1.0f;
  float rot = 0.0f;
  int select = 0;
  if(sscanf(line, "%f %f %d %s %f %f", &t->lambda_min, &t->lambda_step, &t->lambda_num, basename, &b, &rot) < 4)
  {
    fprintf(stderr, "[spectral] could not read hdri file name! usage: <lambda_min> <lambda_step> <lambda_num> <basefilename> <brightness> <rot>\n");
    return 1;
  }
  t->envmap = (envmap_t *)malloc(t->lambda_num * sizeof(envmap_t));
  t->phase = M_PI*rot/180.0;
  t->select = select;
  for(int bin=0;bin<t->lambda_num;bin++)
  {
    snprintf(filename, 1024, "%s_l%d.hdr", basename, (int)(t->lambda_min + t->lambda_step*bin));
    FILE *f = fopen(filename, "rb");
    if(!f)
    {
      char fname[1024];
      snprintf(fname, 1024, "hdri/%s", filename); f = fopen(fname, "rb");
    }
    if(!f)
    {
      fprintf(stderr, "[spectral] could not read hdri file: %s!\n", filename);
      return 1;
    }
    RGBE_ReadHeader(f, &t->width, &t->height, NULL);
    // assume these stay constant over the hdr files:
    if(t->width != 2*t->height) { fprintf(stderr, "[spectral] ERROR: width has to be 2*height!\n"); return 1; }
    t->w2n = 1; t->h2n = 1;
    while(t->w2n<=t->width)  t->w2n <<= 1; t->w2n >>= 1;
    while(t->h2n<=t->height) t->h2n <<= 1; t->h2n >>= 1;
    // printf("wd ht %d %d %d %d %f %f\n", t->width, t->height, t->w2n, t->h2n, t->aspectx, t->aspecty);
    t->aspectx = t->width/(float)t->w2n;
    t->aspecty = t->height/(float)t->h2n;
    t->levels = 0;
    int size = 3*t->width*t->height;
    for(int h = t->h2n;h>0;h>>=1,t->levels++) size += h*h*2;

    t->envmap[bin].pixels = (float *)malloc(sizeof(float)*size);
    RGBE_ReadPixels_RLE(f, t->envmap[bin].pixels, t->width, t->height);
    fclose(f);

    t->envmap[bin].mip = (float **)malloc(sizeof(float *)*(t->levels));
    t->envmap[bin].mip[0] = t->envmap[bin].pixels + 3*t->width*t->height;
    for(int k=0;k<t->levels-1;k++) t->envmap[bin].mip[k+1] = t->envmap[bin].mip[k] + (t->w2n>>k)*(t->h2n>>k);

    // multiply everything by sin theta
    for(int j=0;j<t->height;j++)
    {
      const float sin = sinf(M_PI*(.5f+j)/(float)t->height);
      for(int i=0;i<3*t->width;i++) t->envmap[bin].pixels[i + 3*t->width*j] *= b*sin;
    }

    // create max importance-sample mip-maps with probabilities, smallest 2x1 pixels wide
    t->envmap[bin].sum = 0.0;
    for(int j=0;j<t->h2n;j++)
      for(int i=0;i<t->w2n;i++)
      {
        t->envmap[bin].sum += (t->envmap[bin].mip[0][i + t->w2n*j] = sky_rgbeis_sh(t->envmap[bin].pixels + 3*((int)(i*t->aspectx) + t->width*(int)(j*t->aspecty))));
      }
    t->envmap[bin].sum /= t->w2n*t->h2n;

    for(int m=1;m<t->levels;m++)
    {
      for(int j=0;j<(t->h2n>>m);j++)
      {
        for(int i=0;i<(t->w2n>>m);i++)
        {
          const float all = t->envmap[bin].mip[m-1][(i*2 +     (t->w2n>>(m-1))*j*2)] + t->envmap[bin].mip[m-1][(i*2 +     (t->w2n>>(m-1))*(j*2+1))] +
                            t->envmap[bin].mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*j*2)] + t->envmap[bin].mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*(j*2+1))];
          t->envmap[bin].mip[m][i + (t->w2n>>m)*j] = .25f*all;
          if(all >= 0.001f)
          {
            t->envmap[bin].mip[m-1][(i*2 +     (t->w2n>>(m-1))*j*2    )] /= all;
            t->envmap[bin].mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*j*2    )] /= all;
            t->envmap[bin].mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*(j*2+1))] /= all;
            t->envmap[bin].mip[m-1][(i*2 +     (t->w2n>>(m-1))*(j*2+1))] /= all;
          }
        }
      }
      const float all = t->envmap[bin].mip[t->levels-1][0] + t->envmap[bin].mip[t->levels-1][1];
      t->envmap[bin].mip[t->levels-1][0] /= all;
      t->envmap[bin].mip[t->levels-1][1] /= all;
    }
  }
  printf("[spectral] loaded %d %dx%d hdr files: %s\n", t->lambda_num, t->width, t->height, basename);
  return 0;
}

