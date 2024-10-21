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

#include "corona_common.h"
#include "sampler_common.h"
#include "shader.h"
#include "points.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <dlfcn.h>
#include <strings.h>

rt_t rt; // provide some linkage
__thread rt_tls_t rt_tls;

// provide fake linkage
int  accel_visible(const struct accel_t *b, const ray_t *ray) { return 0; }
void accel_intersect(const struct accel_t *b, const ray_t *ray, hit_t *hit) {}
const float *accel_aabb(const struct accel_t *b) { return 0; }
void accel_print_info(FILE *fd) { }

void writeP5(uint8_t *buf, int size, const char *filename)
{
  FILE *f = fopen(filename, "wb");
  fprintf(f, "P5\n%d %d\n255\n", size, size);
  fwrite(buf, 1, size*size, f);
  fclose(f);
}

void writeP6(uint8_t *buf, int size, const char *filename)
{
  FILE *f = fopen(filename, "wb");
  fprintf(f, "P6\n%d %d\n255\n", size, size);
  fwrite(buf, 1, 3*size*size, f);
  fclose(f);
}

// TODO: really only test for deviation and print images if it doesn't fit.
void write_testimages(shader_so_t *shader, float roughness, int reflect, int count)
{
  const int REFLECT = reflect;
  const int spp = 8;//64;
  int size = 512;
  const int samples = spp * size*size;
  uint8_t buf[size*size];
  float buff[size*size];
  float bufh[size*size];
  char filename[512];
  // int max = 0;
  float maxf = 0.0f, maxh = 0.0f;
  path_t path;
  srand48(666);

  path.index = 0;
  path_init(&path, path.index, 0); // all 0
  path.lambda = 525.0f;
  
  path.length = 2; // now sampling v[2]
  hit_t *hit = &path.v[1].hit;
  hit->n[2] = 1.0f;
  hit->gn[2] = 1.0f;
  hit->prim.extra = 0;
  hit->prim.shapeid = 0;
  hit->prim.vi = 0;
  hit->prim.mb = 0;
  hit->prim.vcnt = 2;
  if(!REFLECT) hit->n[2] = hit->gn[2] = -1.0f;

  path.v[0].hit.prim = INVALID_PRIMID;
  path_volume_vacuum(&path.e[0].vol);
  path.v[0].interior = path.e[0].vol;
  path.e[1].vol = path.e[0].vol;

  vertex_shading_t *sh = &path.v[1].shading;
  sh->rs = 0.06f;
  sh->rd = .8f;
  sh->em = 0.0f;
  sh->rg = 1.0f;
  sh->roughness = roughness;

  get_onb(hit->n, hit->a, hit->b);

  const int num_k = count-1;
  // for(int d=3;d<10;d++)
  // const int d = 3;
  for(int k=0;k<=num_k;k++)
  {
    const float u = k/(float)(num_k+.5f),
                v = 0;//drand48();
    path.e[1].omega[0] = sqrtf(u) * sinf(2.0f*M_PI*v);
    path.e[1].omega[1] = sqrtf(u) * cosf(2.0f*M_PI*v);
    path.e[1].omega[2] = REFLECT ? - sqrtf(1-u) : sqrtf(1-u);

    maxf = 0;
    maxh = 0;
    memset(buff, 0, sizeof(float)*size*size);
    memset(bufh, 0, sizeof(float)*size*size);
    // FIXME: needs decorrelated random samples!
// #pragma omp parallel for firstprivate(path) default(shared)
    for(int i=0;i<samples;i++)
    {
      path.index++;
      shader->prepare(&path, 1, shader->data);
      if(!REFLECT)
        path.e[2].vol = path.v[1].interior;
      else
        path.e[2].vol.ior = 1.0f;

      path.v[1].mode = s_absorb;
      float w = shader->sample(&path, shader->data);
      // TODO: assert omega_out be normalised!
      // normalise(omega_out);
      // sample_cos(omega_out, omega_out+1, omega_out+2, drand48(), drand48());
      // const float w = 1.0f;
      if(!(path.e[2].omega[0] == path.e[2].omega[0]) ||
         !(path.e[2].omega[0] == path.e[2].omega[0]) ||
         !(path.e[2].omega[0] == path.e[2].omega[0]))
      {
        // fprintf(stderr, "sample() outputs NaN!\n");
        path.e[2].omega[0] = path.e[2].omega[0] = 0.0;
        w = 1e20f;
      }
      if(path.e[2].omega[2] <= 0.0f) continue; // wrong hemisphere, not checking that currently.
      if(w <= 0.0f) continue;
      int x = size*.5*(path.e[2].omega[0]+1.);
      int y = size*.5*(path.e[2].omega[1]+1.);
#pragma omp critical
      {
      // buff[x + size*y]++;
      buff[x + size*y] += w/((float)samples);    // estimate bsdf
      bufh[x + size*y] += 1.0f/((float)samples); // track histogram (estimate pdf)
      // buff[x + size*y] *= (1.0f/(float)(size*size))/M_PI;
      maxf = buff[x+size*y] > maxf ? buff[x+size*y] : maxf;
      maxh = bufh[x+size*y] > maxh ? bufh[x+size*y] : maxh;
      }
      // maxf += 1.0/((float)size*size*100);//buff[x + size*y];
    }
    for(int i=0;i<size*size;i++) buf[i] = buff[i]*(255.0f/(float)maxf);
    // sprintf(filename, "estimated_bsdf_%02d_%02d.pgm", d, k);
    sprintf(filename, "%02d_ebsdf.pgm", k+1);
    writeP5(buf, size, filename);

    for(int i=0;i<size*size;i++) buf[i] = bufh[i]*(255.0f/(float)maxh);
    sprintf(filename, "%02d_epdf.pgm", k+1);
    writeP5(buf, size, filename);

    float ebsdf = 0.0f;
    for(int i=0;i<size*size;i++) ebsdf += buff[i];

    float epdf = 0.0f;
    for(int i=0;i<size*size;i++) epdf += bufh[i];

    // commented out: use same scale as histogram.
    // maxf = 0.0f;
    // maxh = 0.0f;
    memset(buff, 0, size*size*sizeof(float));
    memset(bufh, 0, size*size*sizeof(float));
    // bzero(buff, sizeof(float)*size*size);
#pragma omp parallel for firstprivate(path) default(shared) collapse(2)
    for(int j=0;j<size;j++) for(int i=0;i<size;i++)
    {
      path.e[2].omega[0] = 2.0* i/(float)size - 1.;
      path.e[2].omega[1] = 2.0* j/(float)size - 1.;
      const float len2 = path.e[2].omega[0]*path.e[2].omega[0]+path.e[2].omega[1]*path.e[2].omega[1];
      path.e[2].omega[2] = sqrtf(1-len2);
      if(len2 < 1.) for(int k=0;k<spp;k++)
      {
        path.e[1].vol.ior = 1.0f;
        path.index++;
        shader->prepare(&path, 1, shader->data);
        if(!REFLECT)
          path.e[2].vol = path.v[1].interior;
        else
          path.e[2].vol.ior = 1.0f;
        // float w = 1.0f/M_PI; // cos is taken care of by density dictated by i,j loop.
        path.v[1].mode = s_absorb;
        float w = shader->brdf(&path, 1, shader->data);
        w *= 4.0f/(float)(spp*size*size);  // /r^2 since integrating to one over the disk
        buff[i + size*j] += w; // visualize bsdf directly
        // buff[i + size*j] -= w;
        // buff[i + size*j] *= buff[i + size*j];
        // if(buff[i+size*j] > 0.0) printf("buff %f\n", buff[i+size*j]);
// #pragma omp critical
//         {
//         maxf = buff[i + size*j] > maxf ? buff[i + size*j] : maxf;
//         }
        // maxf += buff[i + size*j];

        // XXX FIXME: this is not how it works unfortunately. fresnel from the other side..?
// #define PDF_FLIP
#ifdef PDF_FLIP
        w = shader->pdf(&path, 2, 1, 1, shader->data);
#else
        w = shader->pdf(&path, 1, 1, 2, shader->data);
#endif
        w *= 4.0f/(float)(spp*size*size);  // /r^2 since integrating to one over the disk
        bufh[i + size*j] += w; // visualize pdf directly
// #pragma omp critical
//         {
//         maxh = bufh[i + size*j] > maxh ? bufh[i + size*j] : maxh;
//         }
      }
    }
    for(int i=0;i<size*size;i++) buf[i] = (uint8_t)(buff[i]*(255.0f/maxf));
    sprintf(filename, "%02d_bsdf.pgm", k+1);
    writeP5(buf, size, filename);

    for(int i=0;i<size*size;i++) buf[i] = (uint8_t)(bufh[i]*(255.0f/maxh));
    sprintf(filename, "%02d_pdf.pgm", k+1);
    writeP5(buf, size, filename);

    float bsdf = 0.0f;
    for(int i=0;i<size*size;i++) bsdf += buff[i];

    float pdf = 0.0f;
    for(int i=0;i<size*size;i++) pdf += bufh[i];

    printf("ebsdf-bsdf-epdf-pdf[%d] %f %f %f %f\n", k, ebsdf, bsdf, epdf, pdf);

#if 0 // test pdf
    // commented out: use same scale as histogram.
    maxf = 0.0f;
    for(int j=0;j<size;j++) for(int i=0;i<size;i++)
    {
      path.e[2].omega[0] = 2.0* i/(float)size - 1.;
      path.e[2].omega[1] = 2.0* j/(float)size - 1.;
      const float len2 = path.e[2].omega[0]*path.e[2].omega[0]+path.e[2].omega[1]*path.e[2].omega[1];
      path.e[2].omega[2] = sqrtf(1-len2);
      if(len2 < 1.)
      {
        // float w = 1.0f/M_PI; // cos is taken care of by density dictated by i,j loop.
        float w = shader->pdf(&path, 1, 1, 2, shader->data);
        w *= 4.0f/(float)(size*size);  // /r^2 since integrating to one over the disk
        bufh[i + size*j] -= w;
        bufh[i + size*j] *= bufh[i + size*j];
        // if(buff[i+size*j] > 0.0) printf("buff %f\n", buff[i+size*j]);
        maxf = buff[i + size*j] > maxf ? buff[i + size*j] : maxf;
        // maxf += buff[i + size*j];
      }
    }
    for(int i=0;i<size*size;i++) buf[i] = (uint8_t)(bufh[i]*(255.0f/maxf));
    sprintf(filename, "pdf_diff_%02d_%02d.pgm", d, k);
    writeP5(buf, size, filename);

    maxf = 0.0f;
    for(int i=0;i<size*size;i++) maxf += buff[i];
    printf("pdf %f\n", maxf);
#endif
  }
}

int main(int argc, char *arg[])
{
  rt.points = points_init(1, 1024);
  shader_so_t shader;
  shader.data = NULL;
  float roughness = 0.01f;
  int reflect = 1;
  int count = 6;
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s libtobetested.so [roughness [reflect [count]]]\n", arg[0]);
    exit(1);
  }
  void *handle = dlopen(arg[1], RTLD_LAZY);
  if (!handle)
  {
    fprintf(stderr, "could not open %s!\n", arg[1]);
    exit(1);
  }
  if(argc > 2) roughness = atof(arg[2]);
  if(argc > 3) reflect = atol(arg[3]);
  if(argc > 4) count = atol(arg[4]);
  shader.sample = (sample_t) dlsym(handle, "sample");
  if(dlerror()) goto error;
  shader.brdf = (brdf_t) dlsym(handle, "brdf");
  if(dlerror()) goto error;
  shader.prepare = (prepare_t) dlsym(handle, "prepare");
  if(dlerror()) shader.prepare = NULL;
  shader.cleanup = (cleanup_t) dlsym(handle, "cleanup");
  if(dlerror()) shader.cleanup = NULL;
  shader.pdf = (pdf_t) dlsym(handle, "pdf");
  if(dlerror()) shader.pdf = NULL;
  shader.init = (init_t) dlsym(handle, "init");
  if(dlerror()) shader.init = NULL;
  fprintf(stderr, "[battle-test] init (please type shader definition or pipe to stdin):\n>");
  if(shader.init && shader.init(stdin, &(shader.data))) goto error;
  fprintf(stderr, "thanks.\n");

  write_testimages(&shader, roughness, reflect, count);

  exit(0);
error:
  fprintf(stderr, "could not initialize shader %s!\n", arg[1]);
  exit(255);
}
