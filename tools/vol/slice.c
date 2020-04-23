#include "vol/vol.h"

#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: slice <input.vol>\n");
    exit(1);
  }
  vol_tree_t *tree = vol_open(argv[1]);
  assert(tree);

  // write everything in a series of pfm as slice map:
  for(int l=0;l<VOL_MOTION_SAMPLES;l++)
  {
    char filename[256];
    snprintf(filename, 256, "slices-%02d.pfm", l);
    FILE *f = fopen(filename, "wb");
    fprintf(f, "PF\n%d %d\n-1.0\n", 512, 512);
    for(int jj=0;jj<512;jj++) for(int ii=0;ii<512;ii++)
    {
      int i = ii & 63;
      int j = jj & 63;
      int k = ii / 64 + 8 * (jj / 64);
      float result[3] = {0};
      int empty = vol_sample(tree, i, j, k, 0, s_vol_density, l/(float)VOL_MOTION_SAMPLES, result);
      if(empty) for(int m=0;m<3;m++) result[m] = 0;
      fwrite(result, sizeof(float), 3, f);
    }
    fclose(f);
  }
  exit(0);
}
