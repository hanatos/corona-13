
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main (int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s input.dmap\n", arg[0]);
    exit(1);
  }
  int wd, ht, dp;
  FILE *fdm = fopen(arg[1], "rb");
  int items_read = fscanf(fdm, "DMAP\n%d %d %d", &wd, &ht, &dp);
  if(items_read != 3) goto error;
  fgetc(fdm);
  float *data = (float *)malloc(sizeof(float)*wd*ht*dp);
  if(!data) goto error;
  items_read = fread(data, sizeof(float)*wd*ht*dp, 1, fdm);
  if(items_read != 1) goto error;

  for(int z=0;z<dp;z++)
  {
    char filename[512];
    snprintf(filename, 512, "depth_%04d.pgm", z);
    FILE *f = fopen(filename, "wb");
    fprintf(f, "P5\n%d %d\n255\n", wd, ht);
    for(int k=0;k<wd*ht;k++)
    {
      uint8_t i = data[k + wd*ht*z]/(data[k+wd*ht*z]+1.f)*255;
      fwrite(&i, sizeof(uint8_t), 1, f);
    }
    fclose(f);
  }

error:
  fclose(fdm);
  exit(0);
}
