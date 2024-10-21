#include "framebuffer.h"

int main(int argc, char *argv[])
{
  // rgb2spec_t *r2s = rgb2spec_init("data/xyz2spec.coeff");
  rgb2spec_t *r2s = rgb2spec_init("data/ergb2spec.coeff");
  if(!r2s) exit(1);

  // load pfm file
  FILE *f = fopen(argv[1], "rb");
  if(!f) exit(2);
  int width = 0, height = 0;
  fscanf(f, "PF\n%d %d\n%*[^\n]", &width, &height);
  fgetc(f); // newline
  float *data = malloc(sizeof(float)*3*width*height);
  fread(data, 1, sizeof(float)*3*width*height, f);
  fclose(f);

  int err = fb_tex_from_float(
    r2s, data, width, height, argv[2], 1);

  free(data);
  rgb2spec_cleanup(r2s);
  exit(err);
}
