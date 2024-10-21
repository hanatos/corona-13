#include "framebuffer.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

float srgb_to_linear(float x)
{
  if(x < 0.04045f) return x / 12.92;
  else
  {
    const float a = 0.055;
    return powf((x+a)/(1.f+a), 2.4f);
  }
}

int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    fprintf(stderr, "[ppm2fb] usage: ./ppm2fb input.ppm output.fb\n");
    exit(1);
  }
  rgb2spec_t *r2s = rgb2spec_init("data/ergb2spec.coeff");
  if(!r2s) exit(1);

  // load ppm file
  FILE *f = fopen(argv[1], "rb");
  if(!f) exit(2);
  int width = 0, height = 0;
  fscanf(f, "P6\n%d %d\n255\n", &width, &height);
  float *data = malloc(sizeof(float)*3*width*height);
  uint8_t *data8 = malloc(sizeof(uint8_t)*3*width*height);
  fread(data8, 1, sizeof(uint8_t)*3*width*height, f);
  fclose(f);
  // i totally hate how ppm are always flipped, no matter what you do:
  for(uint64_t j=0;j<height;j++)
    for(uint64_t i=0;i<3*width;i++)
      //data[3*width*(height-1-j)+i] = data8[3*width*j+i]/255.0f;
      data[3*width*(height-1-j)+i] = srgb_to_linear(data8[3*width*j+i]/255.0f);
      // data[3*width*(height-1-j)+i] = powf(data8[3*width*j+i]/255.0f, 1.0f/2.8f);
  free(data8);

  int err = fb_tex_from_float(
    r2s, data, width, height, argv[2], 0);

  free(data);
  rgb2spec_cleanup(r2s);
  exit(err);
}
