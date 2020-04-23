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

#include <cstdio>
#include <cstdlib>
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfRgbaFile.h>

using namespace std;

int main(int argc, char** argv)
{
  char outname[512];
  if (argc < 2)
  {
    printf("Usage: %s infile.pfm [outfile.exr]\n", argv[0]);
    return 1;
  }

  if (argc == 2)
  {
    if (strlen(argv[1]) < 5) {fprintf(stderr, "ERROR: filename too short"); return 0;}
    strncpy(outname, argv[1], sizeof(outname));
    char* extension = outname + strlen(outname) - 4;
    strcpy(extension, ".exr");
  }
  else
  {
    strncpy(outname, argv[2], sizeof(outname));
  }

  FILE *f = fopen(argv[1], "rb");
  if (f)
  {
    int resX, resY;
    float* pData;
    fscanf(f, "PF\n%d %d\n%*[^\n]%*c", &resX, &resY);
    pData = new float[3 * resX * resY];
    fread(pData, 1, sizeof(float)*3*resX*resY, f);
    fclose(f);

    Imf::RgbaOutputFile rgbafile(outname, resX, resY, Imf::WRITE_RGBA);
    Imf::Rgba *pEXRData = new Imf::Rgba[resX * resY];
    for (int j = 0; j < resY; j++)
    for (int i = 0; i < resX; i++)
    {
      pEXRData[i + resX*j].r = fmaxf(0.0f, pData[(i + resX*(j))*3+0]);
      pEXRData[i + resX*j].g = fmaxf(0.0f, pData[(i + resX*(j))*3+1]);
      pEXRData[i + resX*j].b = fmaxf(0.0f, pData[(i + resX*(j))*3+2]);
      pEXRData[i + resX*j].a = 1.0f;
    }
    rgbafile.setFrameBuffer(pEXRData, 1, resX);
    rgbafile.writePixels(resY);

    delete [] pData;
    delete [] pEXRData;
  }
  else
  {
    fprintf(stderr, "ERROR: Could not open %s\n", argv[1]);
    return 0;
  }
  exit(0);
}
