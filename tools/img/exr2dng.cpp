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

using namespace Imf;
using namespace Imath;
using namespace LatLongMap;

// dng code is GPL and stolen :|


#define II       1
#define MM       2  
#define BYTE     1 
#define ASCII    2
#define SHORT    3
#define LONG     4
#define RATIONAL 5
#define SRATIONAL 10

static inline void screenshot_write_buf(uint8_t *buf, int adr, int val)
{
  buf[adr+3] = val & 0xff;
  buf[adr+2] = (val>>8) & 0xff;
  buf[adr+1] = (val>>16) & 0xff;
  buf[adr  ] = val>>24;
}

static inline uint8_t *screenshot_make_tag ( uint16_t tag, uint16_t type, uint32_t lng, 
  uint32_t fld, uint8_t *b, uint8_t *cnt )
{
  screenshot_write_buf(b, 0, (tag<<16)| type);
  screenshot_write_buf(b, 4, lng);
  screenshot_write_buf(b, 8, fld);
  *cnt = *cnt + 1;
  return b + 12;
}

static inline void screenshot_convert_rational(float f, int32_t *num, int32_t *den)
{
  int32_t sign = 1;
  if(f < 0)
  {
    sign = -1;
    f = - f;
  }
  float mult = 1.0f;
  while(f*mult - (int)(f*mult) > 0.0001f) mult++;
  *den = mult;
  *num = (int)(*den * f);
  *num *= sign;
}

static inline void screenshot_write_tiff_header ( FILE *fp, uint32_t xs, uint32_t ys, float Tv, float Av, float f, float iso)
{
  const uint32_t channels = 3;
  uint8_t *b, *offs1, *offs2;
  uint32_t exif_offs;
  uint8_t  buf[1024]; 
  uint8_t  cnt = 0;

  memset(buf, 0, 1024);
  /* TIFF file header.  */
  buf[0] = 0x4d;
  buf[1] = 0x4d;
  buf[3] = 42;
  buf[7] = 10;

  b = buf + 12;
  b = screenshot_make_tag(  254, LONG, 1, 0, b, &cnt ); /* New subfile type.  */
  b = screenshot_make_tag(  256, SHORT, 1, (xs<<16), b, &cnt ); /* Image width.  */
  b = screenshot_make_tag(  257, SHORT, 1, (ys<<16), b, &cnt ); /* Image length.  */
  b = screenshot_make_tag(  258, SHORT, 3, 506, b, &cnt ); /* Bits per sample.  */
  b = screenshot_make_tag(  259, SHORT, 1, (1<<16), b, &cnt ); /* Compression.  */
  b = screenshot_make_tag(  262, SHORT, 1, 34892<<16, b, &cnt);//34892, b, &cnt ); // linear raw /* Photo interp.  */
  b = screenshot_make_tag(  271, ASCII, 8, 494, b, &cnt); // maker, needed for dcraw
  b = screenshot_make_tag(  272, ASCII, 9, 484, b, &cnt); // model
  offs2 = b + 8;
  b = screenshot_make_tag(  273, LONG, 1, 584, b, &cnt ); /* Strip offset.  */
  b = screenshot_make_tag(  277, SHORT, 1, channels<<16, b, &cnt ); /* Samples per pixel.  */
  b = screenshot_make_tag(  278, SHORT, 1, (ys<<16), b, &cnt); /* Rows per strip.  */
  b = screenshot_make_tag(  279, LONG, 1, (ys*xs*channels*2), b, &cnt ); // 16 bits/channel /* Strip byte count.  */
  b = screenshot_make_tag(  284, SHORT, 1, (1<<16), b, &cnt ); /* Planar configuration.  */
  b = screenshot_make_tag(  306, ASCII, 20, 428, b, &cnt ); // DateTime
  offs1 = b + 8;// + 3;
  b = screenshot_make_tag(34665, LONG, 1, 264, b, &cnt); // exif ifd
  b = screenshot_make_tag(50706, BYTE, 4, (1<<24) | (2<<16), b, &cnt); // DNG Version/backward version
  b = screenshot_make_tag(50707, BYTE, 4, (1<<24) | (1<<16), b, &cnt);
  b = screenshot_make_tag(50708, ASCII, 9, 484, b, &cnt); // unique camera model
  b = screenshot_make_tag(50721, SRATIONAL, 9, 328, b, &cnt); // ColorMatrix1 (XYZ->native cam)
  // b = screenshot_make_tag(50728, RATIONAL, 3, 512, b, &cnt); // AsShotNeutral
  b = screenshot_make_tag(50729, RATIONAL, 2, 512, b, &cnt); // AsShotWhiteXY
  b = screenshot_make_tag( 0, 0, 0, 0, b, &cnt ); /* Next IFD.  */
  buf[11] = cnt-1;
  // printf("offset: %d\n", b - buf);
  // set exif IFD offset
  exif_offs = b - buf;
  screenshot_write_buf(buf, offs1 - buf, b - buf);

  b += 2;
  b = screenshot_make_tag(33434, RATIONAL, 1, 400, b, &cnt); // exposure time
  b = screenshot_make_tag(33437, RATIONAL, 1, 408, b, &cnt); // FNumber
  b = screenshot_make_tag(34855, SHORT, 1, ((int)iso)<<16, b, &cnt); // iso speed rating
  b = screenshot_make_tag(37386, RATIONAL, 1, 416, b, &cnt); // focal length
  b = screenshot_make_tag( 0, 0, 0, 0, b, &cnt ); /* Next IFD.  */
  // buf[253] = cnt-buf[11]-1;
  buf[exif_offs + 1] = cnt-buf[11]-1;

  // printf("offset: %d\n", b - buf);

  int32_t num, den;
  // ColorMatrix1 (adobe sRGB D65)
  float m[9] =
  {
      3.24071f,    -0.969258f,    0.0556352f,  
    -1.53726f ,    1.87599f  ,  -0.203996f  , 
    -0.498571f,    0.0415557f,   1.05707f
  };
  for(int k=0;k<9;k++)
  {
    screenshot_convert_rational(m[3*(k%3)+k/3], &num, &den);
    screenshot_write_buf(buf, 328+8*k,   num);
    screenshot_write_buf(buf, 328+8*k+4, den);
  }
  // for(int k=332;k<400;k+=8) screenshot_write_buf(buf, k, 1); // den
  // screenshot_write_buf(buf, 328, 1);// color matrix1: identity
  // screenshot_write_buf(buf, 360, 1);
  // screenshot_write_buf(buf, 392, 1);
  screenshot_convert_rational(Tv, &num, &den);
  screenshot_write_buf(buf, 400, num); // exposure time
  screenshot_write_buf(buf, 404, den);
  screenshot_convert_rational(Av, &num, &den);
  screenshot_write_buf(buf, 408, num); // fnumber
  screenshot_write_buf(buf, 412, den);
  screenshot_convert_rational(f, &num, &den);
  screenshot_write_buf(buf, 416, num); // focal length
  screenshot_write_buf(buf, 420, den);
  strncpy((char *)buf+428, "2008:07:15 13:37:00\0", 20); // DateTime: leet-time
  strncpy((char *)buf+484, "exr2dng0\0", 9);
  strncpy((char *)buf+494, "alfonso\0", 8);

  // bits per sample
  buf[507] = buf[509] = buf[511] = 16; 
  // AsShotNeutral
  // screenshot_convert_rational(0.333, &num, &den);
  // for(int k=0;k<3;k++)
  // {
  //   screenshot_write_buf(buf, 518+8*k,   num);
  //   screenshot_write_buf(buf, 518+8*k+4, den);
  // }
  // AsShotWhiteXY
  screenshot_convert_rational(0.3333, &num, &den);
  screenshot_write_buf(buf, 512, num);
  screenshot_write_buf(buf, 516, den);
  screenshot_convert_rational(0.333, &num, &den);
  screenshot_write_buf(buf, 520, num);
  screenshot_write_buf(buf, 524, den);

  // screenshot_write_buf(buf, offs2-buf, 584);
  fwrite ( buf, 1, 584, fp );
}

static inline void screenshot_write(const char *filename, Imf::Array2D<Imf::Rgba> &pixels, int width, int height, int iso)
{
  uint16_t col;
  FILE* f = fopen(filename, "wb");
  float isof = iso/100.0f;
  if(f)
  {
    screenshot_write_tiff_header(f, width, height, 1./125., 16., 50, 100);
    for(int i=0;i<height;i++)
    {
      for(int k=0;k<width;k++)
      {
        // TODO: use ISO value here?
        col = (uint16_t)(65535*fminf(1.0f, fmaxf(0.0f, isof*pixels[i][k].r)));
        col = (col<<8) | (col >> 8);                               
        fwrite(&col, sizeof(uint16_t), 1, f);                      
        col = (uint16_t)(65535*fminf(1.0f, fmaxf(0.0f, isof*pixels[i][k].g)));
        col = (col<<8) | (col >> 8);                               
        fwrite(&col, sizeof(uint16_t), 1, f);                      
        col = (uint16_t)(65535*fminf(1.0f, fmaxf(0.0f, isof*pixels[i][k].b)));
        col = (col<<8) | (col >> 8);
        fwrite(&col, sizeof(uint16_t), 1, f);
      }
     }
    fclose(f);
  }
}

int main(int argc, char *argv[])
{
  char outname[512];
  if (argc < 3)
  {
    printf("Usage: %s infile.exr iso [outfile.dng]\n", argv[0]);
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
  int iso = atol(argv[2]);
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
  screenshot_write(outname, pixels, width, height, iso);
  exit(0);
}

