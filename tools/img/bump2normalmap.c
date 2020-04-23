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

int main(int argc, char *arg[])
{
  char outname[512];
  float intensity = 1.0f;
  if(argc < 3)
  {
    fprintf(stderr, "usage: %s <intensity> in.pgm [out.ppm]\n", arg[0]);
    exit(1);
  }
  intensity = atof(arg[1]);
  if (argc == 3)
  {
    if (strlen(arg[2]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, arg[2], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".ppm");
  }
  else
  {
    strncpy(outname, arg[3], sizeof(outname));
  }

  FILE *f = fopen(arg[2], "rb");
  if(!f)
  {
    fprintf(stderr, "could not open %s for reading!\n", arg[2]);
    exit(2);
  }
  int width = 0, height = 0;
  fscanf(f, "P5\n%d %d\n%*[^\n]\n", &width, &height);
  unsigned char *idata = (unsigned char *)malloc(sizeof(unsigned char)*width*height);
  unsigned char *odata = (unsigned char *)malloc(sizeof(unsigned char)*3*width*height);
  fread(idata, 1, sizeof(unsigned char)*width*height, f);
  fclose(f);

  f = fopen(outname, "wb");
  if(!f)
  {
    fprintf(stderr, "could not open %s for writing!\n", outname);
    exit(1);
  }

  const int hw = 1;
  const int hh = 1;
  const int du[] = {1,  -2, 1,
                    2,  -4, 2,
                    1,  -2, 1};
  const int dv[] = { 1,  2,  1,
                    -2, -4, -2,
                     1,  2,  1};
  float normal[3];

  const float one255 = 1.0f/255.0f;
  int index = 0;
  for(int j=0;j<height;j++)
  {
    for(int i=0;i<width;i++)
    {
      normal[0] = 0.0f;
      normal[1] = 0.0f;
      if(j>0 && j < height-1 && i>0 && i < width-1) for(int k=-hh;k<=hh;k++) for(int l=-hw;l<=hw;l++)
      {
        normal[0] += intensity * du[k+hw + (2*hw+1)*(l+hh)]*idata[i+k + width*(j+l)]*one255;
        normal[1] += intensity * dv[k+hw + (2*hw+1)*(l+hh)]*idata[i+k + width*(j+l)]*one255;
      }
      normal[2] = 1.0f;
      const float onelen = 1.0f/sqrtf(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      for(int k=0;k<3;k++) normal[k] *= onelen;
      // for(int k=0;k<3;k++) odata[index++] = idata[i + width*j];
      odata[index++] = (unsigned char)fmaxf(0.0f, fminf(255.0f, 255.0f*.5f*(1.0f + normal[0])));
      odata[index++] = (unsigned char)fmaxf(0.0f, fminf(255.0f, 255.0f*.5f*(1.0f + normal[1])));
      odata[index++] = (unsigned char)fmaxf(0.0f, fminf(255.0f, 255.0f*normal[2]));
    }
  }

  fprintf(f, "P6\n%d %d\n255\n", width, height);
  fwrite(odata, 1, 3*width*height, f);
  fclose(f);

  free(idata);
  free(odata);

  exit(0);
}
