#include "camera.h"
#include "corona_common.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  camera_t cam0;
  camera_t cam1;
  if(argc < 3)
  {
    fprintf(stderr, "usage: turntablecam first.cam  inital camera\n"
                    "                    target.cam final frame camera\n"
                    "  [--shutter/-s <degrees>]     shutter degrees, 360 means open up until next frame\n"
                    "  [--frames/-f <n>]            number of frames, camera n == target\n");
    exit(1);
  }
  if(camera_read(&cam0, argv[1]))
  {
    fprintf(stderr, "could not read camera `%s'!\n", argv[1]);
    exit(2);
  }
  if(camera_read(&cam1, argv[2]))
  {
    fprintf(stderr, "could not read camera `%s'!\n", argv[2]);
    exit(2);
  }

  int num_frames = 360;
  float shutter = 180.0f;
  for(int i=3;i<argc;i++)
  {
    if((!strcmp(argv[i], "--shutter") || !strcmp(argv[i], "-s")) && (i+1 < argc))
      shutter = atof(argv[++i]);
    if((!strcmp(argv[i], "--frames") || !strcmp(argv[i], "-f")) && (i+1 < argc))
      num_frames = atol(argv[++i]);
  }

  fprintf(stderr, "writing %d frames with %g degree shutter\n", num_frames, shutter);

  char filename[1024];
  for(int k=0;k<num_frames;k++)
  {
    float t0 = k / (num_frames - 1.0);
    float t1 = (k+(shutter / 360.0)) / (num_frames - 1.0);
    camera_t curr = cam0;
    quaternion_slerp(&cam0.orient, &cam1.orient, t0, &curr.orient);
    curr.orient_t1 = curr.orient; // XXX
    // quaternion_slerp(&cam0.orient_t1, &cam1.orient, t1, &curr.orient);
    for(int i=0;i<3;i++) curr.pos[i]    = (1.0-t0) * cam0.pos[i] + t0 * cam1.pos[i];
    for(int i=0;i<3;i++) curr.pos_t1[i] = (1.0-t1) * cam0.pos[i] + t1 * cam1.pos[i];
    curr.focus     = (1.0-t0) * cam0.focus     + t0 * cam1.focus;
    curr.focal_length = (1.0-t0) * cam0.focal_length + t0 * cam1.focal_length;
    curr.iso = (1.0-t0) * cam0.iso + t0 * cam1.iso;
    curr.aperture_value = (1.0-t0) * cam0.aperture_value + t0 * cam1.aperture_value;
    curr.exposure_value = (1.0-t0) * cam0.exposure_value + t0 * cam1.exposure_value;
    snprintf(filename, 1024, "lerpcam_%04d.cam", k);
    if(camera_write(&curr, filename))
    {
      fprintf(stderr, "could not write camera `%s'\n", filename);
      exit(4);
    }
  }
  exit(0);
}
