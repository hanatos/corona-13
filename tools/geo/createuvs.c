#include "corona_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

float *read_input(const char *filename, size_t *cnt)
{
  FILE *fd = fopen(filename, "rb");

  if(!fd) return NULL;
  fseek (fd, 0, SEEK_END);
  const size_t filesize = ftell(fd);
  fseek (fd, 0, SEEK_SET);
  float *uv = (float *)malloc(filesize);
  int items = fread(uv, filesize, 1, fd);
  if(items != 1) exit(2);
  fclose(fd);
  *cnt = filesize;
  return uv;
}

int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    fprintf(stderr, "usage: %s input.uv input.ra2\n", argv[0]);
    exit(1);
  }

  size_t cnt;
  float *all_uv = read_input(argv[1], &cnt);
  const int num = cnt/(2*sizeof(float));
  float *all_v = read_input(argv[2], &cnt);
  if(num != cnt/(3*sizeof(float)))
  {
    fprintf(stderr, "uv and ra2 don't contain the same number of points!\n");
    exit(3);
  }

  // get bounding box
  float aabb[6] = {FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};
  for(int k=0;k<num;k++)
  {
    for(int i=0;i<3;i++) aabb[i] = fminf(aabb[i], all_v[3*k+i]);
    for(int i=0;i<3;i++) aabb[3+i] = fmaxf(aabb[i+3], all_v[3*k+i]);
  }

  // optimized for material tester probe, to fall into the triangular cutout
  // TODO: make that seam go away in abovementioned cutout and where phi wraps!
  // float pole[3] = {-1.0, 1.0, -2.0};
  float pole[3] = {.0, -1.0, .0};
  const float scale = 10.0f;
  normalise(pole);
  float a[3], b[3];
  get_onb(pole, a, b);
  float *uv = all_uv, *v = all_v;
  for(int prim=0;prim<num/3;prim++)
  {
    float trin[3];
    trin[0] = (v[3*1+1] - v[3*0+1])*(v[3*2+2] - v[3*0+2]) - (v[3*1+2] - v[3*0+2])*(v[3*2+1] - v[3*0+1]);
    trin[1] = (v[3*1+2] - v[3*0+2])*(v[3*2+0] - v[3*0+0]) - (v[3*1+0] - v[3*0+0])*(v[3*2+2] - v[3*0+2]);
    trin[2] = (v[3*1+0] - v[3*0+0])*(v[3*2+1] - v[3*0+1]) - (v[3*1+1] - v[3*0+1])*(v[3*2+0] - v[3*0+0]);
    // if(uv[0] == 0.0f && uv[1] == 0.0f)
    for(int vtx=0;vtx<3;vtx++)
    {
      // get polar coordinates based on the vertices
      float d[3] = {0.0f};
      for(int i=0;i<3;i++) d[i] = v[i] - (aabb[i] + aabb[3+i])*.5f;
      normalise(d);

      const float uu = dotproduct(d, a);
      const float vv = dotproduct(d, b);
      const float ww = dotproduct(d, pole);
      const float phi = atan2f(uu, vv);
      const float theta = acosf(fminf(1.0f, fmaxf(-1.0f, ww)));
      if(dotproduct(d, trin) > 0.0f)
        uv[0] = -scale * phi/M_PI;
      else
        uv[0] = scale * phi/M_PI;
      uv[1] = scale * theta/M_PI;
      // fprintf(stderr, "uv %f %f\n", uv[0], uv[1]);
      uv += 2;
      v += 3;
    }
  }
  // write new uvs:
  FILE *fd = fopen("output.uv", "wb");
  cnt = fwrite(all_uv, num*2*sizeof(float), 1, fd);
  fclose(fd);

  free(all_uv);
  free(all_v);

  exit(0);
}
