#pragma once

#include "quaternion.h"
#include "corona_common.h"
#include "mf.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>

// struct camera_t;
// typedef struct camera_t camera_t;
// camera description, also straight file format when written to disk:
typedef struct camera_t
{
  char magic[4];             // magic number, has to be the four letters 'CCAM'
  int32_t version;           // version number
  float pos[3];              // world space position
  float pos_t1[3];           // shutter close version
  quaternion_t orient;       // orientation of frame
  quaternion_t orient_t1;    // shutter close version

  float speed;               // world space movement speed for gui

  float focus_sensor_offset; // offset of sensor from focus-at-infinity position, in [mm]
  float focus;               // focus distance, in [dm]=1.0f in world space

  float film_width;          // film back width in [mm]
  float film_height;         // film back height in [mm]
  float crop_factor;         // scale factor of film back wrt full frame (35mm)
  int aperture_value;        // aperture, access view_f_stop[.] to get f/# number
  int exposure_value;        // exposure, access view_exposure_time[.] to get [s]
  float focal_length;        // focal length, if applicable
  float iso;                 // iso value of the film back.
}
camera_t;

struct path_t;
typedef struct path_t path_t;

// output the camera info on the display
void camera_display_info(const camera_t *c);

// output camera info for render sidecar file
void camera_print_info(const camera_t *c, FILE *f);

// set camera's focus distance (in decimeter, 1.0f units)
void camera_set_focus(camera_t *c, float dist);

// sample initial point on the given path
mf_t camera_sample(const camera_t *c, path_t *p);

// evaluate pdf of sampling vertex v (0 or p->length-1)
mf_t camera_pdf(const camera_t *c, const path_t *p, int v);

// connect path to camera, either updating vertex p->v[0] or the last one,
// depending on tracing direction (p->v[0].mode & s_sensor?)
mf_t camera_connect(const camera_t *c, path_t *p);

// pdf of connecting given vertex to given camera as p->v[v]
mf_t camera_pdf_connect(const camera_t *c, const path_t *p, int v);

// mutate point on aperture
mf_t camera_mutate_aperture(const camera_t *c, path_t *p, const float r1, const float r2, const float step);

// pdf of mutating point on aperture curr->tent
mf_t camera_pdf_mutate_aperture(const camera_t *c, const path_t *curr, const path_t *tent, const float step);

// evaluate camera responsivity for given path
mf_t camera_eval(const camera_t *c, const path_t *p);

// global initialisation
int camera_global_init();

// file io

// support loading of legacy cameras
typedef struct camera_v0_t
{
  int32_t legacy0;
  float pos[3];
  quaternion_t orient;
  float speed;
  int32_t legacy1[7];
  float iso;
  quaternion_t orient_t1; // orientation and position
  float pos_t1[3];        // for motion blur
  float focus_sensor_offset;
  float fill[4];
  float focus;
  float legacy2;
  float crop_factor;
  float film_width;
  float film_height;
  int aperture_value;
  float focal_length;
  float legacy3;
  int exposure_value;
}
camera_v0_t;

static inline int camera_write(const camera_t *cam, const char *const filename)
{
  FILE* f = fopen(filename, "wbx");
  if(f)
  {
    fwrite(cam, sizeof(camera_t), 1, f);
    fclose(f);
    return 0;
  }
  return 1;
}

static inline int camera_write_humanly_readable(const camera_t *cam, const char *const filename)
{
  FILE* f = fopen(filename, "wbx");
  if(f)
  {
    hit_t hit;
    hit.b[1] = hit.n[2] = 1.0f;
    hit.b[0] = hit.b[2] = hit.n[0] = hit.n[1] = 0.0f;
    const quaternion_t *q = &(cam->orient);
    quaternion_transform(q, hit.b);
    quaternion_transform(q, hit.n);
    normalise(hit.b);
    normalise(hit.n);

    fprintf(f, "pos_world = (%f, %f, %f)\n",cam->pos[0], cam->pos[1], cam->pos[2]);
    fprintf(f, "lookat = (%f, %f, %f)\n", cam->pos[0] + hit.n[0], cam->pos[1] + hit.n[1], cam->pos[2] + hit.n[2]);
    fprintf(f, "up = (%f, %f, %f)\n", hit.b[0], hit.b[1], hit.b[2]);

    fprintf(f, "frame_orientation w = %f, X = (%f, %f, %f)\n", cam->orient.w, cam->orient.x[0], cam->orient.x[1], cam->orient.x[2]);

    fprintf(f, "focus_sensor_offset = %f\n", cam->focus_sensor_offset);
    fprintf(f, "focus = %f\n", cam->focus);

    fprintf(f, "film width = %f\n", cam->film_width);
    fprintf(f, "film height = %f\n", cam->film_height);
    fprintf(f, "crop_factor = %f\n", cam->crop_factor);
    fprintf(f, "aperture_value = %d\n", cam->aperture_value);
    fprintf(f, "exposure_value = %d\n", cam->exposure_value);
    fprintf(f, "focal_length = %f\n", cam->focal_length);
    //field of view in degree for x and y axis, can be used by mitsuba
    fprintf(f, "fov_x = %f\n", atan2(cam->film_width * 0.5f, cam->focal_length) * 360.f / M_PI);
    fprintf(f, "fov_y = %f\n", atan2(cam->film_height * 0.5f, cam->focal_length) * 360.f / M_PI);
    fprintf(f, "iso = %f\n", cam->iso);

    fclose(f);
    return 0;
  }
  return 1;
}

static inline int camera_read(camera_t *cam, const char *const filename)
{
  FILE* f = fopen(filename, "rb");
  if(f)
  {
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    int bytes_read = 0;
    fseek(f, 0, SEEK_SET);
    if(size == sizeof(camera_v0_t))
    { // load legacy camera
      camera_v0_t c0;
      bytes_read = fread(&c0, sizeof(camera_v0_t), 1, f);
      if(bytes_read != 1)
      {
        fclose(f);
        return 1;
      }
      for(int k=0;k<3;k++) cam->pos[k] = c0.pos[k];
      for(int k=0;k<3;k++) cam->pos_t1[k] = c0.pos_t1[k];
      cam->orient = c0.orient;
      cam->orient_t1 = c0.orient_t1;
      cam->speed = c0.speed;
      cam->focus = c0.focus;
      cam->film_width  = c0.film_width;
      cam->film_height = c0.film_height;
      cam->aperture_value = c0.aperture_value;
      cam->exposure_value = c0.exposure_value;
      cam->focal_length = c0.focal_length;
      cam->iso = c0.iso;
    }
    else if(size == sizeof(camera_t))
    {
      bytes_read = fread(cam, sizeof(camera_t), 1, f);
      if(bytes_read != 1 || strncmp(cam->magic, "CCAM", 4) || cam->version != 1)
      {
        fclose(f);
        return 1;
      }
    }
    return 0;
  }
  return 1;
}
