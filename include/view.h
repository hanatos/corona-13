#pragma once
#include "corona_common.h"
#include "quaternion.h"
#include "mf.h"
#include <stdint.h>

// declare our structs of unknown content
struct view_t;
typedef struct view_t view_t;

// fwd declare
struct path_t;
typedef struct path_t path_t;

// fwd declare
struct camera_t;
typedef struct camera_t camera_t;

// move active camera
typedef enum view_move_t
{
  s_view_move_fw = 1<<0,
  s_view_move_bk = 1<<1,
  s_view_move_rg = 1<<2,
  s_view_move_lf = 1<<3,
  s_view_move_up = 1<<4,
  s_view_move_dn = 1<<5,
}
view_move_t;

// change parameters of active camera
typedef enum view_ctl_t
{
  s_view_ctl_av_up,      // aperture value
  s_view_ctl_av_dn,
  s_view_ctl_fl_up,      // focal length
  s_view_ctl_fl_dn,
  s_view_ctl_tv_up,      // time value
  s_view_ctl_tv_dn,
  s_view_ctl_speed_up,   // movement speed
  s_view_ctl_speed_dn,
  s_view_ctl_iso_up,     // iso
  s_view_ctl_iso_dn,
}
view_ctl_t;

view_t *view_init();
void view_cleanup(view_t *v);


// frame buffer and path space interfaces:
// =======================================

// reset all frame buffers and sampler and pointsampler
void view_clear();

// reset only the frame buffers
void view_clear_frame(view_t *r);

// sample camera id from [0,1) random number
int view_sample_camid(const float rand);

// pdf to sample this camera id
float view_pdf_camid(const int camid);

// splat given path with spectral contribution
void view_splat(const struct path_t *path, const mf_t value);

// splat given path with spectral contribution and arbitrary gaussian
void view_splat_gaussian(const struct path_t *path, const float mu[2], const float S[4], const float col[3]);

// don't splat, but convert to colour and accumulate path length stats.
void view_deferred_splat(const struct path_t *path, const mf_t value, float *col);

// splat given path with camera space colour col
void view_splat_col(const struct path_t *path, const float *col);

// output render to file
void view_write_images(const char *filename);

// render a progression, colormanage, and push result to display
void view_render();

// push colourmanaged pixels to display
void view_update();

// sample time that seems appropriate for the shutter we have.
// need to call view_sample_camid on the path before this.
float view_sample_time(const struct path_t *p, const float rand);

// pdf to have sampled time given the shutter of the camera used on this path
float view_pdf_time(const struct path_t *p, const float time);



// wrap camera interface functions:
// ================================

// output the camera info on the display
void view_cam_display_info();

// output descriptive info to sidecar file
void view_print_info(FILE *fd);

// sample initial point on the given path
mf_t view_cam_sample(struct path_t *p);

// evaluate pdf of sampling vertex v (0 or p->length-1)
mf_t view_cam_pdf(const struct path_t *p, int v);

// connect path to camera, either updating vertex p->v[0] or the last one,
// depending on tracing direction (p->v[0].mode & s_sensor?)
mf_t view_cam_connect(struct path_t *p);

// pdf of connecting given vertex to given camera as p->v[v]
mf_t view_cam_pdf_connect(const struct path_t *p, int v);

// mutate point on aperture
mf_t view_cam_mutate_aperture(struct path_t *p, const float r1, const float r2, const float step);

// pdf of mutating point on aperture curr->tent
mf_t view_cam_pdf_mutate_aperture(const struct path_t *curr, const struct path_t *tent, const float step);

// evaluate camera responsivity for given path
mf_t view_cam_eval(const struct path_t *p);

// initialise frame for camera and time associated with path
void view_cam_init_frame(const path_t *p, hit_t *hit);

// helper to convert time value to time in seconds
float view_tv2time(int tv);

// helper to convert aperture value to f-stop
float view_av2fstop(int av);

// access frame buffer dimensions
uint64_t view_width();
uint64_t view_height();

// return number of progressions rendered so far
uint64_t view_overlays();

// set exposure gain. this is useful for metropolis methods
// which compute a histogram that needs to be scaled by mean image brightness.
void view_set_exposure_gain(const float gain);

// get camera
const camera_t *view_get_camera(const path_t *path);

// gui related functions:
// ======================

void view_set_active_camid(int camid);
int view_get_active_camid();

// gui interaction, movement:
// ==========================

// rotate active camera
void view_active_cam_rotate(const float *axis, const float angle);

// apply movement indicated by move_begin and move_end
void view_active_cam_move();

// move camera (key down)
void view_move_begin(view_move_t m);

// move camera (key up)
void view_move_end(view_move_t m);

// change camera parameters
void view_ctl(view_ctl_t c);

// focus active camera to object under given pixel.
// for stereo, focuses both left and right eye and adjust convergence.
void view_set_focus(float i, float j, int verbose);

// print cmdline options
void view_print_usage();

// write currently active camera to disk
int view_cam_write(const char *const filename);

// write currently active camera to disk
// in a format that can be read by humans
int view_cam_write_readable(const char *const filename_readable);

// read currently active camera from file
int view_cam_read(const char *const filename);
