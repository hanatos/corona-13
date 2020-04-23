#include "camera.h"
#include "corona_common.h"
#include <stdio.h>
#include <stdlib.h>

void init_frame(const camera_t *c, const float time, hit_t *hit)
{
  quaternion_t q;
  quaternion_slerp(&c->orient, &c->orient_t1, time, &q);
  hit->a[0] = hit->b[1] = hit->n[2] = 1.0f;
  hit->a[1] = hit->a[2] = hit->b[0] = hit->b[2] = hit->n[0] = hit->n[1] = 0.0f;
  quaternion_transform(&q, hit->a);
  quaternion_transform(&q, hit->b);
  quaternion_transform(&q, hit->n);
  for(int k=0;k<3;k++) hit->gn[k] = hit->n[k];
  for(int k=0;k<3;k++) hit->x[k] = c->pos[k] * (1.0f-time) + c->pos_t1[k] * time;
  normalise(hit->a);
  normalise(hit->b);
  normalise(hit->n);
}

int main(int argc, char *argv[])
{
  camera_t cam;
  if(argc < 4)
  {
    fprintf(stderr, "usage: turntablecam first.cam  inital camera\n"
                    "  [--shutter/-s <degrees>]     shutter degrees, 360 means open up until next frame\n"
                    "  [--frames/-f <n>]            number of frames for full circle, camera n == camera 0\n");
    exit(1);
  }
  if(camera_read(&cam, argv[1]))
  {
    fprintf(stderr, "could not read camera `%s'!\n", argv[1]);
    exit(2);
  }

  int num_frames = 360;
  float shutter = 180.0f;
  for(int i=2;i<argc;i++)
  {
    if((!strcmp(argv[i], "--shutter") || !strcmp(argv[i], "-s")) && (i+1 < argc))
      shutter = atof(argv[++i]);
    if((!strcmp(argv[i], "--frames") || !strcmp(argv[i], "-f")) && (i+1 < argc))
      num_frames = atol(argv[++i]);
  }

  hit_t hit;
  init_frame(&cam, 0.0f, &hit);
  const float dist = cam.focus;
  const float pivot[3] = {
    cam.pos[0] + hit.n[0]*dist,
    cam.pos[1] + hit.n[1]*dist,
    cam.pos[2] + hit.n[2]*dist};

  float up[3] = {0.0, 0.0, 1.0};
  quaternion_t quat = cam.orient;
  quaternion_invert(&quat);
  quaternion_transform(&quat, up); // worldspace up in camera space for rotation

  char filename[1024];
  for(int k=0;k<num_frames;k++)
  {
    camera_t curr = cam;
    const float angle0 = (360.0f*k)/num_frames * M_PI/180.0f;
    const float angle1 = (360.0f*(k+(shutter/360.0)))/num_frames * M_PI/180.0f;
    quaternion_t rot;
    quaternion_init(&rot, angle0, up);
    quaternion_mult(&(curr.orient), &rot);
    quaternion_init(&rot, angle1, up);
    quaternion_mult(&(curr.orient_t1), &rot);
    init_frame(&curr, 0.0f, &hit);
    for(int i=0;i<3;i++) curr.pos[i] = pivot[i] - hit.n[i]*dist;
    init_frame(&curr, 1.0f, &hit);
    for(int i=0;i<3;i++) curr.pos_t1[i] = pivot[i] - hit.n[i]*dist;
    snprintf(filename, 1024, "turntable_%04d.cam", k);
    if(camera_write(&curr, filename))
    {
      fprintf(stderr, "could not write camera `%s'\n", filename);
      exit(4);
    }
  }
  exit(0);
}
