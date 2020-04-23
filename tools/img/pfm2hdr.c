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


int main(int argc, char *arg[])
{
  char outname[512];
  if(argc < 2)
  {
    fprintf(stderr, "%s: usage in.pfm [out.hdr]\n", arg[0]);
    exit(1);
  }
  if (argc == 2)
  {
    if (strlen(arg[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, arg[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".hdr");
  }
  else
  {
    strncpy(outname, arg[2], sizeof(outname));
  }

  FILE *f = fopen(arg[1], "rb");
  if(!f)
  {
    fprintf(stderr, "could not open %s for reading!\n", arg[1]);
    exit(2);
  }
  int width = 0, height = 0;
  fscanf(f, "PF\n%d %d\n%*[^\n]\n", &width, &height);
  float *data = (float *)malloc(3*sizeof(float)*width*height);
  fread(data, 1, 3*sizeof(float)*width*height, f);
  for(int k=0;k<3*width*height;k++) data[k] = fmaxf(0.0f, data[k]);
  fclose(f);

  f = fopen(outname, "wb");
  if(!f)
  {
    fprintf(stderr, "could not open %s for reading!\n", arg[2]);
    exit(2);
  }
  RGBE_WriteHeader(f, width, height, NULL);
  RGBE_WritePixels(f, data, width*height);
  fclose(f);
  free(data);
  exit(0);
}

