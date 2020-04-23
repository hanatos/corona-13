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


int main(int argc, char *argv[])
{
  char outname[512];
  if (argc < 3)
  {
    printf("usage: %s infile.pfm subtract.pfm [outfile.pfm]\n", argv[0]);
    return 1;
  }

  if (argc == 3)
  {
    if (strlen(argv[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, argv[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, "_diff.pfm");
  }
  else
  {
    strncpy(outname, argv[3], sizeof(outname));
  }
  uint32_t width, height, wd, ht;

  FILE *fin = fopen(argv[1], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[1]); exit(1); }
  fscanf(fin, "PF\n%d %d\n%*[^\n]", &width, &height);
  fgetc(fin); // \n
  float *pixels = (float *)malloc(width*height*3*sizeof(float));
  fread(pixels, width*height*3, sizeof(float), fin);
  fclose(fin);

  fin = fopen(argv[2], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[2]); exit(2); }
  fscanf(fin, "PF\n%d %d\n%*[^\n]", &wd, &ht);
  if(wd != width || ht != height) { fprintf(stderr, "image dimensions do not match! %dx%d vs %dx%d\n", width, height, wd, ht); exit(3); }
  fgetc(fin); // \n
  float *diff = (float *)malloc(width*height*3*sizeof(float));
  fread(diff, width*height*3, sizeof(float), fin);
  fclose(fin);

  FILE *f = fopen(outname, "wb");
  if(f)
  {
    char header[1024];
    snprintf(header, 1024, "PF\n%d %d\n-1.0", width, height);
    size_t len = strlen(header);
    fprintf(f, "PF\n%d %d\n-1.0", width, height);
    size_t off = 0;
    while((len + 1 + off) & 0xf) off++;
    while(off-- > 0) fprintf(f, "0");
    fprintf(f, "\n");
    double rmse = 0.0f;
    for(uint64_t k=0;k<(uint64_t)width*height;k++)
    {
      float c[3] = {
        fabsf(pixels[3*k]-diff[3*k]),
        fabsf(pixels[3*k+1]-diff[3*k+1]),
        fabsf(pixels[3*k+2]-diff[3*k+2])};
      rmse += c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
      fwrite(&c, sizeof(float), 3, f);
    }
    rmse = sqrt(rmse/((double)width*height));
    fprintf(stdout, "[pfmdiff] rmse: %g\n", rmse);
  }
  else
  {
    fprintf(stderr, "could not open %s for writing!\n", outname);
    exit(4);
  }

  free(pixels);
  free(diff);
  exit(0);
}

