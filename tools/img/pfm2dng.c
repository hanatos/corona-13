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
#include "../include/screenshot_dng.h"


int main(int argc, char *argv[])
{
  char outname[512];
  if (argc < 3)
  {
    printf("compiled with\n");
    colorspace_print_info(stdout);
    printf("Usage: %s infile.pfm iso [outfile.dng]\n", argv[0]);
    return 1;
  }

  if (argc == 3)
  {
    if (strlen(argv[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, argv[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".dng");
  }
  else
  {
    strncpy(outname, argv[3], sizeof(outname));
  }
  float iso = atof(argv[2]);
  float *pixels, drggn;
  int width;
  int height;

  FILE *fin = fopen(argv[1], "rb");
  if(!fin) { fprintf(stderr, "could not open %s!\n", argv[1]); exit(1); }
  fscanf(fin, "PF\n%d %d\n%*[^\n]%*c", &width, &height);
  // fscanf(fin, "PF\n%d %d\n%f\n", &width, &height, &drggn);
  // fseek(fin, 17, SEEK_SET);
  // printf("w, h, order %d %d %f\n", width, height, drggn);
  // printf("offset: %d\n", ftell(fin));
  // fscanf(fin, "PF\n%d %d\n-1.0\n", &width, &height);
  pixels = (float *)malloc(width*height*3*sizeof(float));
  fread(pixels, width*height*3, sizeof(float), fin);
#if 0
  for(int k=0;k<width*height*3;k++)
  {
    int dreggn = *(int *)(pixels + k);
    dreggn = (dreggn << 24) | ((dreggn & 0xff00) << 8) | ((dreggn & 0xff0000) >> 8) | (dreggn >> 24);  
    pixels[k] = *(float *)&dreggn;
  }
#endif
  fclose(fin);

  // FILE *dreggn = fopen("dreggnaoeu.ppm", "wb");
  // fprintf(dreggn, "P6\n%d %d\n255\n", width, height);
  // for(int k=0;k<width*height*3;k++) { char c = (int)(pixels[k]/(pixels[k]+1.0f)); fwrite(&c, 1, 1, dreggn); }
  // fclose(dreggn);

  uint16_t col[3];
  FILE* f = fopen(outname, "wb");
  if(f)
  {
    screenshot_write_tiff_header(f, width, height, 0.0f, 0.0f, 0.0f, iso);
    iso /= 100.0f;
    for(int i=0;i<height;i++)
    {
      for(int k=0;k<width;k++)
      {
        int p = 3*(width*i + k);
        for(int j=0;j<3;j++) { col[j] = (uint16_t)(65535*fminf(1.0f, fmaxf(0.0f, iso*pixels[p+j]))); col[j] = (col[j]<<8) | (col[j] >> 8); }
        fwrite(col, sizeof(uint16_t), 3, f);
      }
     }
    fclose(f);
  }
  else
  {
    fprintf(stderr, "could not open %s for writing!\n", outname);
    exit(2);
  }

  free(pixels);
  exit(0);
}

