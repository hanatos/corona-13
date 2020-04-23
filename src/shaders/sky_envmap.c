/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "shader.h"
#include "accel.h"
#include "framebuffer.h"
#include "sampler_common.h"

typedef struct
{
  float *pixels;
  float **mip;
  double sum;
  int levels;
  int width, w2n;
  int height, h2n;
  float aspectx, aspecty;
  float world[9];
  float world_inv[9];
  framebuffer_t fb;
  float mul; // brightness scale
}
rgbe_t;

float sky_envmap_sh(const float *const coeff)
{
  return coeff[3];// XXX
  const __m128 l4 = _mm_set_ps(400.0f, 480.0f, 560.0f, 660.0f);
  const __m128 eval = _mm_mul_ps(rgb2spec_eval_sse(coeff, l4), _mm_set1_ps(coeff[3]));
  return eval[0]+eval[1]+eval[2]+eval[3];
}

void cleanup(void *data)
{
  rgbe_t *t = data;
  free(t->pixels);
  free(t->mip);
  fb_cleanup(&t->fb);
  free(t);
}

// use y up:
// #define FLIP // was 2 0 "" for anim
#define FLIP_A 1
#define FLIP_B 2
#define FLIP_SIGN -1.0*

mf_t eval(path_t *p, int v, void *data)
{
  rgbe_t *t = data;
  float dir_w[3];
  if(v == 0)
    for(int k=0;k<3;k++) dir_w[k] = -p->e[1].omega[k];
  else
    for(int k=0;k<3;k++) dir_w[k] = p->e[v].omega[k];
#ifdef FLIP
  // flip backward
  const float tmp = dir_w[FLIP_A];
  dir_w[FLIP_A] = FLIP_SIGN dir_w[FLIP_B];
  dir_w[FLIP_B] = tmp;
#endif
  float dir[3] = { 0 };
  mat3_mulv(t->world_inv, dir_w, dir);
  float x, y;
  if(fabsf(dir[2]) > 1.0f)
  {
    y = 0;
    x = 0;
  }
  else
  {
    y = acosf(dir[2])/M_PI * t->height;
    x = ((M_PI + atan2f(dir[0], dir[1]))/(2*M_PI)) * t->width;
  }
  if(!(x < t->width && x >= 0.0 && y < t->height && y >= 0.0)) x = y = 0.0;
  const int i = (int)x, j = (int)y;
  const float *cf = fb_fetchi(&t->fb, i, j);
  return mf_mul(mf_rgb2spec(cf, p->lambda), mf_set1(t->mul * cf[3]));
}

mf_t sample(path_t *p, void *data)
{
  rgbe_t *t = (rgbe_t *)data;
  float x1, x2;
  if(p->v[0].mode & s_emit)
  { // light ray
    x1 = pointsampler(p, s_dim_edf_x);
    x2 = pointsampler(p, s_dim_edf_y);
  }
  else
  { // next event
    x1 = pointsampler(p, s_dim_nee_x);
    x2 = pointsampler(p, s_dim_nee_y);
  }
  double wx = x1, wy = x2;

  if(wx < t->mip[t->levels-1][0])
    wx *= .5/t->mip[t->levels-1][0];
  else
    wx = 1.0 - (1.0-wx)*.5/(1.0 - t->mip[t->levels-1][0]);
  wx *= 2.0;
  for(int m=t->levels-2;m>=0;m--)
  {
    int i = (int)wx, j = (int)wy, wd = t->w2n>>m;
    if(2*i >= wd) i = (wd-1)/2;
    if(4*j >= wd) j = (wd-1)/4;
    if(wy < j + (t->mip[m][2*i + 2*j*wd]+t->mip[m][2*i+1 + 2*j*wd]))
    {
      wy = j + (wy-j)*.5/(t->mip[m][2*i + 2*j*wd]+t->mip[m][2*i+1 + 2*j*wd]);
      j = 2*j;
    }
    else if (t->mip[m][2*i + 2*j*wd]+t->mip[m][2*i+1 + 2*j*wd] < 0.999)
    {
      wy = j + 1.0 - (1.0-(wy-j))*.5/(1.0 - (t->mip[m][2*i + 2*j*wd]+t->mip[m][2*i+1 + 2*j*wd]));
      j = 2*j + 1;
    }
    const double f = t->mip[m][2*i + j*wd]/(t->mip[m][2*i + j*wd] + t->mip[m][2*i+1 + j*wd]);
    if(wx <= i + f)
      wx = i + (wx-i)*.5/f;
    else if(wx > i + f) // do nothing for nan
      wx = i + 1.0 - (1.0-wx+i)*.5/(1.0-f);
    wx *= 2.0; wy *= 2.0f;
  }

  const float x = wx*t->aspectx, y = wy*t->aspecty;
  const int i = (int)x, j = (int)y;
  const float *coeff = fb_fetchi(&t->fb, i, j);
  const float theta = M_PI*(y/t->height);
  const float phi   = 2*M_PI*(x/t->width) - M_PI;

  float sin_theta, cos_theta;
  sincosf(theta, &sin_theta, &cos_theta);
  const float quantized_sin_theta = sinf(M_PI*(.5f+j)/t->height);

  float sin_phi, cos_phi;
  sincosf(phi, &sin_phi, &cos_phi);
  const int v = p->length;
  p->e[v].omega[0] = sin_phi*sin_theta;
  p->e[v].omega[1] = cos_phi*sin_theta;
  p->e[v].omega[2] = cos_theta;
  p->e[v].dist = FLT_MAX;
  p->v[v].shading.roughness = 1.0f;
  p->v[v].flags = s_environment;
  p->v[v].mode = s_emit;
  p->v[v].pdf = mf_set1(sky_envmap_sh(coeff)*t->mul*quantized_sin_theta/(M_PI*M_PI*2.0*t->sum*sin_theta));
  p->v[v].shading.em = mf_mul(mf_rgb2spec(coeff, p->lambda),
      mf_set1(t->mul * coeff[3]));

  float dir_tmp[3] = { p->e[v].omega[0], p->e[v].omega[1], p->e[v].omega[2] };
  mat3_mulv(t->world, dir_tmp, p->e[v].omega);
  float *dir = p->e[v].omega;
#ifdef FLIP
  // flip forward
  const float tmp = dir[FLIP_A];
  dir[FLIP_A] = dir[FLIP_B];
  dir[FLIP_B] = FLIP_SIGN tmp;
#endif
  if(p->length)
  {
    const float *aabb = accel_aabb(rt.accel);
    const float far = aabb[3] + aabb[4] + aabb[5]
      - aabb[0] - aabb[1] - aabb[2];
    for(int k=0;k<3;k++)
      p->v[v].hit.x[k] = p->v[v-1].hit.x[k] + far*dir[k];
  }
  for(int k=0;k<3;k++)
    p->v[v].hit.n[k] = p->v[v].hit.gn[k] = - dir[k];
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = -1;
  return mf_div(p->v[v].shading.em, p->v[v].pdf);
}

mf_t pdf(path_t *p, int v, void *data)
{
  rgbe_t *t = data;
  float dir_w[3];
  if(p->v[0].mode & s_emit)
    for(int k=0;k<3;k++) dir_w[k] = -p->e[v].omega[k];
  else
    for(int k=0;k<3;k++) dir_w[k] = p->e[v].omega[k];
#ifdef FLIP
  // flip backward
  const float tmp = dir_w[FLIP_A];
  dir_w[FLIP_A] = FLIP_SIGN dir_w[FLIP_B];
  dir_w[FLIP_B] = tmp;
#endif
  float dir[3] = { 0 };
  mat3_mulv(t->world_inv, dir_w, dir);
  const float y = CLAMP(acosf(dir[2])/M_PI * t->height, 0, t->height-1);
  const float x = CLAMP(((M_PI + atan2f(dir[0], dir[1]))/(2*M_PI)) * t->width, 0, t->width-1);
  const float sin_theta = sqrtf(MAX(1e-12f, 1.0f-dir[2]*dir[2]));
  if(!(x < t->width && x >= 0.0 && y < t->height && y >= 0.0)) assert(0);//x = y = 0.0;
  const int i = (int)x, j = (int)y;
  const float quantized_sin_theta = sinf(M_PI*(.5f+j)/t->height);
  return mf_set1(sky_envmap_sh(fb_fetchi(&t->fb, i, j))*t->mul*quantized_sin_theta/(t->sum * sin_theta * 2.0f*M_PI*M_PI));
}

#if 0
float test_pdf( void *data)
{
  srand48(666);
  path_t p;
  path_init(&p, 0, 0);
  p.lambda = 550.f;
  int N = 1<<20;
  double sum = 0.0;

  for(int i=0;i<N;i++)
  {
    const float r0 = drand48(), r1 = drand48();
    sample_sphere(p.e[1].omega, p.e[1].omega+1, p.e[1].omega+2, r0, r1);
    const double f = pdf(&p, 1, data) * 4.0f * M_PI;
    // fprintf(stderr, "rand %g wo %g pdf = %g\n", r0,
        // p.e[0].omega[0], // p.e[0].omega[1],
        // p.e[0].omega[2], f);
    sum += f;
  }
  fprintf(stderr, "mean pdf = %g\n", sum/N);

  p.length = 2;
  sum = 0.0;
  for(int i=0;i<N;i++)
  {
    const float r0 = drand48(), r1 = drand48();
    sample_sphere(p.e[1].omega, p.e[1].omega+1, p.e[1].omega+2, r0, r1);
    const double f = eval(&p, data) * 4.0f * M_PI;
    // fprintf(stderr, "rand %g wo %g pdf = %g\n", r0,
        // p.e[0].omega[0], // p.e[0].omega[1],
        // p.e[0].omega[2], f);
    sum += f;
  }
  fprintf(stderr, "mean eval (spherical) = %g\n", sum/N);

  sum = 0.0;
  p.length = 1;
  for(int i=0;i<N;i++)
  {
    // const float r0 = drand48(), r1 = drand48();
    const double f = sample(&p, data);
    const double g = pdf(&p, 1, data);
    // assert(fabs(p.v[1].pdf - g)<1e-3f);
    sum += f;
  }
  fprintf(stderr, "mean eval (importance sample) = %g\n", sum/N);
  exit(0);
}
#endif

int init(FILE* s, void **data)
{
  rgbe_t *t = malloc(sizeof(rgbe_t));
  *data = t;
  char filename[1024];
  char line[1024];
  int dreggn = 0;
  dreggn = fscanf(s, "%[^\n]\n", line);
  float b = 1.0f;

  float rot_x = 0.f, rot_y = 0.f, rot_z = 0.f;
  if(sscanf(line, "%s %f %f %f %f", filename, &b, &rot_x, &rot_y, &rot_z) < 1)
  {
	  if(sscanf(line, "%s %f", filename, &b) < 1)
	  {
		fprintf(stderr, "[envmap] could not read hdri file name! usage: <filename> <brightness> [rot_x] [rot_y] [rot_z]\n");
		return 1;
	  }
  }
  t->mul = b;

  float rx[9] = {0};
  float ry[9] = {0};
  float rz[9] = {0};
  float rtemp[9] = {0};
  float ax[3] = {1, 0, 0};
  float ay[3] = {0, 1, 0};
  float az[3] = {0, 0, 1};
  mat3_rotate(ax, rot_x, rx);
  mat3_rotate(ay, rot_y, ry);
  mat3_rotate(az, rot_z, rz);
  mat3_mul(ry, rz, rtemp);
  mat3_mul(rx, rtemp, t->world);
  mat3_invert(t->world, t->world_inv);

  if(fb_map(&t->fb, filename))
  {
    fprintf(stderr, "[envmap] could not read hdri file: %s!\n", filename);
    return 1;
  }

  t->width  = t->fb.header->width;
  t->height = t->fb.header->height;

  if(t->width != 2*t->height) { fprintf(stderr, "[envmap] ERROR: width has to be 2*height!\n"); return 1; }
  t->w2n = 1; t->h2n = 1;
  while(t->w2n<=t->width)  t->w2n <<= 1;
  t->w2n >>= 1;
  while(t->h2n<=t->height) t->h2n <<= 1;
  t->h2n >>= 1;
  t->aspectx = t->width/(float)t->w2n;
  t->aspecty = t->height/(float)t->h2n;
  // printf("[envmap] aspect %g %g\n", t->aspectx, t->aspecty);
  // if(_popcnt32(t->width) != 1) { fprintf(stderr, "[envmap] ERROR: width has to be a power of 2!\n"); return 1; }
  // if(_popcnt32(t->height) != 1) { fprintf(stderr, "[envmap] ERROR: height has to be a power of 2!\n"); return 1; }
  t->levels = 0;
  int size = 0;
  for(int h = t->h2n;h>0;h>>=1,t->levels++) size += h*h*2;
  t->pixels = common_alloc(16, sizeof(float)*size);

  t->mip = (float **)malloc(sizeof(float *)*(t->levels));
  t->mip[0] = t->pixels;
  for(int k=0;k<t->levels-1;k++) t->mip[k+1] = t->mip[k] + (t->w2n>>k)*(t->h2n>>k);

  // create max importance-sample mip-maps with probabilities, smallest 2x1 pixels wide
  t->sum = 0.0;
  for(int j=0;j<t->h2n;j++) for(int i=0;i<t->w2n;i++)
  {
    const int ii = (i * t->aspectx + 0.5f);
    const int jj = (j * t->aspecty + 0.5f);
    const float sin = sinf(M_PI*jj/(float)t->height);
    t->mip[0][i + t->w2n*j] = sky_envmap_sh(fb_fetchi(&t->fb, ii, jj)) * b*sin;
    t->sum += t->mip[0][i+t->w2n*j];
  }

  t->sum /= t->w2n*t->h2n;

  for(int m=1;m<t->levels;m++)
  {
    for(int j=0;j<(t->h2n>>m);j++)
    {
      for(int i=0;i<(t->w2n>>m);i++)
      {
        const float all = t->mip[m-1][(i*2 +     (t->w2n>>(m-1))*j*2)] + t->mip[m-1][(i*2 +     (t->w2n>>(m-1))*(j*2+1))] +
                          t->mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*j*2)] + t->mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*(j*2+1))];
        t->mip[m][i + (t->w2n>>m)*j] = .25f*all;
        if(all >= 0.001f)
        {
          t->mip[m-1][(i*2 +     (t->w2n>>(m-1))*j*2    )] /= all;
          t->mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*j*2    )] /= all;
          t->mip[m-1][(i*2 + 1 + (t->w2n>>(m-1))*(j*2+1))] /= all;
          t->mip[m-1][(i*2 +     (t->w2n>>(m-1))*(j*2+1))] /= all;
        }
      }
    }
    const float all = t->mip[t->levels-1][0] + t->mip[t->levels-1][1];
    t->mip[t->levels-1][0] /= all;
    t->mip[t->levels-1][1] /= all;
  }
  printf("[envmap] loaded %dx%d hdr file: %s\n", t->width, t->height, filename);
  return dreggn == -1;
}

