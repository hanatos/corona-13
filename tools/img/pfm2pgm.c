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
#include <string.h>


int main(int argc, char *argv[])
{
  char outname[512];
  if (argc < 3)
  {
    printf("Usage: %s infile.pfm iso [outfile.pgm]\n", argv[0]);
    return 1;
  }

  if (argc == 3)
  {
    if (strlen(argv[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, argv[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".pgm");
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
  // fscanf(fin, "PF\n%d %d\n%*[^\n]%*c", &width, &height);
  fscanf(fin, "%d %d\n", &width, &height);
  // fscanf(fin, "PF\n%d %d\n%f\n", &width, &height, &drggn);
  // fseek(fin, 17, SEEK_SET);
  // printf("w, h, order %d %d %f\n", width, height, drggn);
  // printf("offset: %d\n", ftell(fin));
  // fscanf(fin, "PF\n%d %d\n-1.0\n", &width, &height);
  pixels = (float *)malloc(width*height*3*sizeof(float));
  for(int k=0;k<width*height;k++)
    fscanf(fin, "%f\n", pixels+k);
  // fread(pixels, width*height*3, sizeof(float), fin);
  fclose(fin);

  FILE *f = fopen(outname, "wb");
  if(f)
  {
    fprintf(f, "P5\n%d %d\n255\n", width, height);
    for(int k=0;k<width*height;k++)
    {
      unsigned char c = pixels[k] * 255;
      fwrite(&c, sizeof(char), 1, f);
    }
  }
  else
  {
    fprintf(stderr, "could not open %s for writing!\n", outname);
    exit(2);
  }

  free(pixels);
  exit(0);
}

