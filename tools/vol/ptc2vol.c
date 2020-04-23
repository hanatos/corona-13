// this little tool converts blender smoke point cache files to hierarchical grid .vol files.

#include "vol/vol.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

static uint32_t *res = 0;
static float *density[VOL_MOTION_SAMPLES] = {0};
static float *temperature[VOL_MOTION_SAMPLES] = {0};
static float gaabb[6] = {0};

// callback used by vol interface to sample the tree
int payload_fill(void *data, vol_payload_uncompressed_t *payload, const float *aabb)
{ // struct comes in memset(0), only write what's needed.
  int have_data = 0;
  for(int k=0;k<512;k++)
  {
    vol_index_t vx = {0};
    vx.idx = k;
    // sample input, fill leaf struct
    float x[3];
    x[0] = aabb[0] + (aabb[0+3]-aabb[0])*(vx.i+.5f)/8.0f;
    x[1] = aabb[1] + (aabb[1+3]-aabb[1])*(vx.j+.5f)/8.0f;
    x[2] = aabb[2] + (aabb[2+3]-aabb[2])*(vx.k+.5f)/8.0f;

    int i[3];
    for(int m=0;m<3;m++) i[m] = (int)((x[m] - gaabb[m])/(gaabb[3+m]-gaabb[m])*(res[m]*res[3]));
    for(int m=0;m<3;m++) if(i[m] < 0 || i[m] >= res[m]*res[3]) continue; // out of bounds, empty
    // compute voxel index, flip z axis (seems to be upside down)
    const int index = res[0]*res[3]*(res[1]*res[3]*i[2] + i[1]) + i[0];

    for(int s=0;s<VOL_MOTION_SAMPLES;s++)
    {
      const float temp = temperature[s][index];
      float dens = density[s][index];
      // increase our chances of empty voxels, simulation has a lot of 1e-30 and lower:
      if(dens < 1e-3f) dens = 0.0f;

      if(dens > 0.0) have_data = 1;
      assert(dens == dens && dens < FLT_MAX && dens >= 0);
      payload->d[k][s] = dens;
      // convert meaningless blender ignition and max temp values to kelvin:
      payload->t[k][s] = 1000*(1.25 + (1.75 - 1.25) * temp);
    }
  }
  return have_data;
}

int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "[ptc2vol] convert blender full grid dumps (.ptc) to hierarchical grids (.vol)\n");
    fprintf(stderr, "[ptc2vol] usage: %s input.ptc%%04d\n", arg[0]);
    exit(1);
  }

  char filename[2048];
  int beg = 9999, end = 0;
  for(int i=0;i<=9999;i++)
  {
    snprintf(filename, sizeof(filename), arg[1], i);
    struct stat sb;
    memset(&sb, 0, sizeof(sb));
    stat(filename, &sb);
    if(S_ISREG(sb.st_mode))
    {
      if(beg == 9999) beg = i;
      end = i;
    }
    else if(beg < 9999) break;
  }
  fprintf(stderr, "[ptc2vol] resampling frames %d to %d\n", beg, end);
  assert(end >= beg + VOL_MOTION_SAMPLES-1);

  int fd[VOL_MOTION_SAMPLES];
  size_t data_size = 0;
  float *data[VOL_MOTION_SAMPLES];
  float voxel_size, *loc, *rot;
  size_t len_hires;

  for(int t=0;t<VOL_MOTION_SAMPLES;t++)
  {
    snprintf(filename, sizeof(filename), arg[1], beg+t);
    fd[t] = open(filename, O_RDONLY);
    if(fd[t] < 0)
    {
      fprintf(stderr, "[ptc2vol] could not open file `%s'!\n", filename);
      exit(2);
    }
    if(!t)
      data_size = lseek(fd[t], 0, SEEK_END);
    else assert(data_size == lseek(fd[t], 0, SEEK_END));
    lseek(fd[t], 0, SEEK_SET);
    data[t] = (float *)mmap(0, data_size, PROT_READ, MAP_SHARED, fd[t], 0);
    if(t)
    { // check that res is still the same!
      for(int k=0;k<4;k++) assert(res[k] == ((uint32_t *)data[t])[k]);
    }
    else
    {
      res = ((uint32_t *)data[0]);
      len_hires = res[0]*res[1]*res[2]*res[3]*res[3]*res[3];
      fprintf(stderr, "[ptc2vol] reading a %dx%dx%d grid, hires * %d\n", res[0], res[1], res[2], res[3]);
      voxel_size = data[0][4];
      loc = data[0] + 5;
      rot = data[0] + 8;
    }
    density[t] = data[t] + 11;
    temperature[t] = data[t] + 11 + len_hires;
  }

  // TODO: get content_box by stepping through the voxels once in advance!

  char output[512];
  snprintf(output, 512, "%s", arg[1]);
  char *c = output + strlen(output);
  for(;c > output && *c != '.';c--);
  snprintf(c, 5, ".vol");

  voxel_size /= res[3];
  fprintf(stderr, "[ptc2vol] transform: voxel size %g off %g %g %g rot %g %g %g\n", voxel_size, loc[0], loc[1], loc[2], rot[0], rot[1], rot[2]);
  fprintf(stderr, "[ptc2vol] creating `%s'\n", output);

  for(int k=0;k<3;k++) gaabb[k] = 0;
  for(int k=3;k<6;k++) gaabb[k] = voxel_size * res[k-3]*res[3];

  if(vol_create_tree(output, gaabb, voxel_size, loc, rot, &payload_fill, 0)) exit(1);

  for(int t=0;t<VOL_MOTION_SAMPLES;t++)
  {
    munmap(data[t], data_size);
    close(fd[t]);
  }

  exit(0);
}
