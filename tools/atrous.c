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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rgbe.h"

void write_pfm(const char *outname, const int w, const int h, const float *buf)
{
  FILE *f = fopen(outname, "wb");
  int ret = 0;
  if(!f)
  {
    fprintf(stderr, "could not open %s for writing\n", outname);
    exit(2);
  }
  fprintf(f, "PF\n%d %d\n-1.0\n", w, h);
  ret = fwrite(buf, 3*sizeof(float), w*h, f);
  fclose(f);
}

#define IND3(i, j) (3*width*(j) + 3*(i))
#define MIN(a,b) ((a)>(b)?(b):(a))

float edge(const float *a, const float *b)
{
  return 1.0/(0.01 + (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1]) + (b[2] - a[2])*(b[2] - a[2]));
}

int main(int argc, char *arg[])
{
  char outname[512];
  if(argc < 2)
  {
    fprintf(stderr, "%s: usage in.pfm\n", arg[0]);
    exit(1);
  }
  if (argc == 2)
  {
    if (strlen(arg[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, arg[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    *extension = '\0';
  }

  FILE *f = fopen(arg[1], "rb");
  if(!f)
  {
    fprintf(stderr, "could not open %s for reading!\n", arg[1]);
    exit(2);
  }
  char filename[1024];
  int width = 0, height = 0, ret = 0;
  ret = fscanf(f, "PF\n%d %d\n%*[^\n]\n", &width, &height);
  float *data = (float *)malloc(3*sizeof(float)*width*height);
  float *tmp  = (float *)malloc(3*sizeof(float)*width*height);
  ret = fread(data, 1, 3*sizeof(float)*width*height, f);
  for(int k=0;k<3*width*height;k++) data[k] = fmaxf(0.0f, data[k]);
  fclose(f);

  float *edges = (float *)malloc(3*sizeof(float)*width*height);
  for(int j=0;j<height-1;j++) for(int i=0;i<width-1;i++)
  {
    edges[IND3(i,j) + 0] = edge(data + IND3(i,j), data + IND3(i+1,j));
    edges[IND3(i,j) + 1] = edge(data + IND3(i,j), data + IND3(i,j+1));
    edges[IND3(i,j) + 2] = edge(data + IND3(i,j), data + IND3(i+1,j+1));
  }

  // copy borders
  const int rad = 3;
  memcpy(tmp, data, sizeof(float)*3*width*height);
  for(int d=0;d<5;d++)
  {
    // blur
    for(int j=rad;j<height-rad;j++) for(int i=rad;i<width-rad;i++)
    {
      for(int k=0;k<3;k++) tmp[IND3(i,j)+k] = 0.0f;
      float weight = 0.0f;
      for(int l=-rad;l<=rad;l++) for(int m=-rad;m<=rad;m++)
      {
        float w = edge(data + IND3(i,j), data + IND3(i+l,j+m));
        //(l == 0 && m == 0) ? w = 1.0 : edges[IND3(MIN(i, i+l), MIN(j,j+m))+((l!=0?1:0)|(m!=0?2:0))-1];
        for(int k=0;k<3;k++) tmp[IND3(i,j)+k] += w*data[IND3(i+l,j+m)+k];
        weight += w;
      }
      for(int k=0;k<3;k++) tmp[IND3(i,j)+k] *= (1.0/weight);
    }
    // write blurry version
    snprintf(filename, 1024, "%s_level_%d_coarse.pfm", outname, d);
    write_pfm(filename, width, height, tmp);
    //recover details
    float max[3] = {0., 0., 0.};
    for(int j=0;j<height;j++) for(int i=0;i<width;i++) for(int k=0;k<3;k++)
    {
      data[IND3(i,j)+k] -= tmp[IND3(i,j)+k];
      data[IND3(i,j)+k] = fabsf(data[IND3(i,j)+k]);
      max[k] = fmaxf(max[k], data[IND3(i,j)+k]);
    }
    for(int j=0;j<height;j++) for(int i=0;i<width;i++) for(int k=0;k<3;k++) data[IND3(i,j)+k] /= max[k];
    snprintf(filename, 1024, "%s_level_%d_fine.pfm", outname, d);
    write_pfm(filename, width, height, data);
    // prepare next iteration.
    memcpy(data, tmp, sizeof(float)*3*width*height);
  }

  free(edges);
  free(data);
  free(tmp);
  exit(0);
}

