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

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void write_pfm(const char *filename, float *buf, uint64_t width, uint64_t height)
{
  FILE *f = fopen(filename, "wb");
  if(f)
  {
    char header[1024];
    snprintf(header, 1024, "PF\n%lu %lu\n-1.0", width, height);
    size_t len = strlen(header);
    fprintf(f, "PF\n%lu %lu\n-1.0", width, height);
    size_t off = 0;
    while((len + 1 + off) & 0xf) off++;
    while(off-- > 0) fprintf(f, "0");
    fprintf(f, "\n");
    fwrite(buf, width*height*3, sizeof(float), f);
    fclose(f);
  }
}

int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    printf("usage: %s primal.pfm gradx.pfm grady.pfm]\n", argv[0]);
    return 1;
  }

  uint64_t width, height, wd, ht;

  FILE *fin = fopen(argv[1], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[1]); exit(1); }
  fscanf(fin, "PF\n%lu %lu\n%*[^\n]", &width, &height);
  fgetc(fin); // \n
  float *pixels = (float *)malloc(width*height*3*sizeof(float));
  fread(pixels, width*height*3, sizeof(float), fin);
  fclose(fin);

  fin = fopen(argv[2], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[2]); exit(2); }
  fscanf(fin, "PF\n%lu %lu\n%*[^\n]", &wd, &ht);
  if(wd != width || ht != height) { fprintf(stderr, "image dimensions do not match! %lux%lu vs %lux%lu\n", width, height, wd, ht); exit(3); }
  fgetc(fin); // \n
  float *gradx = (float *)malloc(width*height*3*sizeof(float));
  fread(gradx, width*height*3, sizeof(float), fin);
  fclose(fin);

  fin = fopen(argv[3], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[3]); exit(4); }
  fscanf(fin, "PF\n%lu %lu\n%*[^\n]", &wd, &ht);
  if(wd != width || ht != height) { fprintf(stderr, "image dimensions do not match! %lux%lu vs %lux%lu\n", width, height, wd, ht); exit(4); }
  fgetc(fin); // \n
  float *grady = (float *)malloc(width*height*3*sizeof(float));
  fread(grady, width*height*3, sizeof(float), fin);
  fclose(fin);

// #pragma omp parallel for schedule(static)
  for(uint64_t j=0;j<height;j++)
    for(uint64_t i=0;i<width;i++)
      for(uint64_t c=0;c<3;c++)
  {
    uint64_t k = 3*(j*width+i)+c;
    const float primal = (pixels[k] + gradx[k] + grady[k])/3.0f;
    if(!(pixels[k] > 0.0)) pixels[k] = 0;
    if(i < width-3)  gradx[k] = gradx[k+3] - pixels[k];
    else             gradx[k] = 0.0;
    if(j < height-3) grady[k] = grady[k+3*width] - pixels[k];
    else             grady[k] = 0.0;
    pixels[k] = primal;
    if(!(pixels[k] > 0.0)) pixels[k] = 0;
  }

  write_pfm("img.pfm", pixels, width, height);
  write_pfm("img_grad_x.pfm", gradx, width, height);
  write_pfm("img_grad_y.pfm", grady, width, height);

  free(pixels);
  free(gradx);
  free(grady);
  exit(0);
}

