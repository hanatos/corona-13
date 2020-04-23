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

#include "shader.h"
#include "rgbe.h"
#include "spectrum.h"

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

typedef struct
{
  float sundir[3];
  float *pixels;
  int width;
  int height;
}
rgbe_t;

float eval(path_t *path, void* p)
{
  float col[3], dir[3];

  if(path->v[0].mode & s_emit)
    for(int k=0;k<3;k++) dir[k] = -path->e[1].omega[k];
  else
    for(int k=0;k<3;k++) dir[k] = path->e[path->length-1].omega[k];
  rgbe_t *t = (rgbe_t *)p;
  float x, y;
  if(fabsf(dir[2]) > 0.999f)
  {
    y = 0;
    x = 0;
  }
  else
  {
    y = acosf(dir[2])/M_PI * (t->height-1);
    x = ((M_PI + atan2f(dir[0], dir[1]))/(2*M_PI)) * (t->width-1);
  }
  if(!(x < t->width && x >= 0.0 && y < t->height && y >= 0.0)) x = y = 0.0;
  const int i = (int)x, j = (int)y;
  /*if(!(x < t->width)) printf("x = %f dir = %f %f %f\n", x, dir[0], dir[1], dir[2]);
  assert(x < t->width);
  if(!(y < t->height)) printf("y = %f dir = %f %f %f\n", y, dir[0], dir[1], dir[2]);
  assert(y < t->height);
  if(!(x >= 0)) printf("x = %f dir = %f %f %f\n", x, dir[0], dir[1], dir[2]);
  assert(x >= 0);
  if(!(y >= 0)) printf("y = %f dir = %f %f %f\n", y, dir[0], dir[1], dir[2]);
  assert(y >= 0);*/
#if 0
  col[0] = t->pixels[3*(t->width*j+i)];
  col[1] = t->pixels[3*(t->width*j+i)+1];
  col[2] = t->pixels[3*(t->width*j+i)+2];
#else
  const float u = x - i, v = y - j;
  const int jp = j == t->height - 1 ? j : j+1;
  const int ip = i == t->width  - 1 ? i : i+1;
  col[0] = (1-u)*((1-v)*t->pixels[3*(t->width*j+i)] + v*t->pixels[3*(t->width*jp+i)]) +
      u*((1-v)*t->pixels[3*(t->width*j+ip)] + v*t->pixels[3*(t->width*jp+ip)]);
  col[1] = (1-u)*((1-v)*t->pixels[3*(t->width*j+i)+1] + v*t->pixels[3*(t->width*jp+i)+1]) +
      u*((1-v)*t->pixels[3*(t->width*j+ip)+1] + v*t->pixels[3*(t->width*jp+ip)+1]);
  col[2] = (1-u)*((1-v)*t->pixels[3*(t->width*j+i)+2] + v*t->pixels[3*(t->width*jp+i)+2]) +
      u*((1-v)*t->pixels[3*(t->width*j+ip)+2] + v*t->pixels[3*(t->width*jp+ip)+2]);

  return spectrum_rgb_to_p(path->lambda, col);
#endif
}

int init(FILE* s, void **data)
{
  rgbe_t *t = (rgbe_t *)malloc(sizeof(rgbe_t));
  *data = t;
  char filename[1024];
  char line[1024];
  int dreggn = 0;
  dreggn = fscanf(s, "%[^\n]\n", line);
  float b = 1.0f;
  for(int k=0;k<3;k++) t->sundir[k] = INFINITY;
  if(sscanf(line, "%s %f %f %f %f", filename, &b, t->sundir, t->sundir+1, t->sundir+2) < 1)
  {
    fprintf(stderr, "[rgbe::init] could not read hdri file name! usage: <filename> <brightness> [<sundir[3]>]\n");
    return 1;
  }
  FILE *f = fopen(filename, "rb");
  if(!f)
  {
    char fname[1024];
    snprintf(fname, 256, "hdri/%s", filename); f = fopen(fname, "rb");
  }
  if(!f)
  {
    fprintf(stderr, "[rgbe::init] could not read hdri file: %s!\n", filename);
    return 1;
  }
  RGBE_ReadHeader(f, &t->width, &t->height, NULL);
  t->pixels = (float *)malloc(sizeof(float)*3*t->width*t->height);
  RGBE_ReadPixels_RLE(f, t->pixels, t->width, t->height);
  fclose(f);
  printf("[rgbe::init] loaded %dx%d hdr file: %s\n", t->width, t->height, filename);
  for(int i=0;i<t->height*t->width*3;i++) t->pixels[i] *= b;

#if 0
  f = fopen("dreggn.ppm", "wb");
  FILE *f2 = fopen("dreggn2.ppm", "wb");
  fprintf(f, "P6\n%d %d\n255\n", t->width, t->height);
  fprintf(f2, "P6\n%d %d\n255\n", t->width, t->height);
  float h[3] = {0,0,0}, h2[3] = {0,0,0};
  const int cnt = 100;
  for(int i=0;i<t->height*t->width*3;i+=3)
  {
    float rgb[3] = {0, 0, 0};
    float rgb2[3] = {0, 0, 0};
    unsigned char buf[3] = {0, 0, 0};
    for(int k=0;k<cnt;k++)
    {
      const float lambda = spectrum_sample_lambda(rand() / (RAND_MAX + 1.0));
      const float power = spectrum_rgb_to_p(lambda, t->pixels+i);
      spectrum_p_to_camera(lambda, power, rgb);
      for(int j=0;j<3;j++) rgb2[j] += rgb[j];
    }
    for(int j=0;j<3;j++) rgb2[j] /= (float)cnt;
    //printf("brightness: %f %f %f\n", rgb2[0], rgb2[1], rgb2[2]);
    for(int j=0;j<3;j++) h[j] += rgb2[j];
    for(int j=0;j<3;j++) h2[j] += t->pixels[i+j];
    for(int j=0;j<3;j++) buf[j] = (unsigned char)(fmaxf(0.0, fminf(1.0f, rgb2[j]/(rgb2[j]+1.0)))*255.0f);
    fwrite(buf, 3, sizeof(unsigned char), f);
    for(int j=0;j<3;j++) buf[j] = (unsigned char)(fmaxf(0.0, fminf(1.0f, t->pixels[i+j]/(t->pixels[i+j]+1.0)))*255.0f);
    fwrite(buf, 3, sizeof(unsigned char), f2);
  }
  for(int j=0;j<3;j++)
  {
    h[j] /= (t->height*t->width);
    h2[j] /= (t->height*t->width);
  }
  printf("brightness: %f %f %f real: %f %f %f\n", h[0], h[1], h[2], h2[0], h2[1], h2[2]);
  fclose(f);
  fclose(f2);
#endif
  return dreggn == -1;
}

#if 0
float sky_daylight_add_sun(const rayhit_t *hit, const ray_t *ray, void *data)
{
  shader_daylight_t* s = (shader_daylight_t *)data;
  const float gamma = acosf(fminf(-dotproduct(ray->dir, s->sun_dir), 1.0f));
  if (gamma < (M_PI*6.312f/(2.0f*180.0f)))
  {
    int i = hit->lambda/10.0f - 38.0f;
    return s->sun_XYZ[1]*(shader_S0_380_780[i] + s->sun_M1*shader_S1_380_780[i] + s->sun_M2*shader_S2_380_780[i]);
  }
  return 0.0f;
}

float sky_daylight_sample_sun(rayhit_t *hit, ray_t *ray, const float x1, const float x2, void *data)
{
  shader_daylight_t* s = (shader_daylight_t *)data;
  float u[3], v[3];
  get_perpendicular(s->sun_dir, u);
  normalise(u);
  crossproduct(s->sun_dir, u, v);
  normalise(v);
  // sample solid angle of sun (approximate as disc)
  const float max_gamma = (M_PI*6.312f/(2.0f*180.0f));
  const float far = rt.accel->aabb[3] + rt.accel->aabb[4] + rt.accel->aabb[5]
    - rt.accel->aabb[0] - rt.accel->aabb[1]- rt.accel->aabb[2];
  const float r = sqrtf(x1)*sinf(max_gamma);
  const float inv_p = M_PI*sinf(max_gamma)*sinf(max_gamma);
  for(int k=0;k<3;k++) ray->dir[k] = -cosf(max_gamma)*s->sun_dir[k] + sinf(2*M_PI*x2)*r*u[k] + cosf(2*M_PI*x2)*r*v[k];
  for(int k=0;k<3;k++) ray->dir[k] *= far;
  int i = hit->lambda/10.0f - 38.0f;
  const float power = s->sun_XYZ[1]*(shader_S0_380_780[i] + s->sun_M1*shader_S1_380_780[i] + s->sun_M2*shader_S2_380_780[i]);
  return inv_p*power;
}
#endif
