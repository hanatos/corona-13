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
#include <ImfEnvmap.h>
#include <half.h>
#include <ImfRgba.h>
#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <ImathBox.h>
#include <inttypes.h>
#include "rgbe.c"

using namespace Imf;
using namespace Imath;
using namespace LatLongMap;

static inline void screenshot_write(const char *filename, Imf::Array2D<Imf::Rgba> &pixels, int width, int height)
{
  uint16_t col;
  FILE* f = fopen(filename, "wb");
  if(f)
  {
    RGBE_WriteHeader(f, width, height, NULL);
    float *data = (float *)malloc(sizeof(float)*3*width*height);
    for(int i=0;i<height;i++)
    {
      for(int k=0;k<width;k++)
      {
        data[3*width*i + 3*k  ] = pixels[i][k].r;
        data[3*width*i + 3*k+1] = pixels[i][k].g;
        data[3*width*i + 3*k+2] = pixels[i][k].b;
      }
     }
    RGBE_WritePixels(f, data, width*height);
    free(data);
    fclose(f);
  }
}

int main(int argc, char *argv[])
{
  char outname[512];
  if (argc < 2)
  {
    printf("Usage: %s infile.exr [outfile.hdr]\n", argv[0]);
    return 1;
  }

  if (argc == 2)
  {
    if (strlen(argv[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, argv[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".hdr");
  }
  else
  {
    strncpy(outname, argv[2], sizeof(outname));
  }
  Imf::Array2D<Imf::Rgba> pixels;
  int width;
  int height;
  RgbaInputFile file (argv[1]);
  Box2i dw = file.dataWindow();
  width = dw.max.x - dw.min.x + 1;
  height = dw.max.y - dw.min.y + 1;
  pixels.resizeErase (height, width);
  file.setFrameBuffer (&(pixels)[0][0] - dw.min.x - dw.min.y * width, 1, width);
  file.readPixels (dw.min.y, dw.max.y);
  screenshot_write(outname, pixels, width, height);
  exit(0);
}

