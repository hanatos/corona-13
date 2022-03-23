#include "view.h"
#include "camera.h"
#include "lights.h"
#include "pathspace.h"
#include "pointsampler.h"
#include "render.h"
#include "sampler.h"
#include "threads.h"
#define FILTER_ATOMIC   // needed for spline filter even in the absence of lt :(
#include "filter.h"
#include "filter/gaussian.h"
#include "filter/box.h"
#include "filter/blackmanharris.h"
#include "framebuffer.h"
#include "screenshot.h"
#include "spectrum.h"
#include "prims.h"

#include <xmmintrin.h>
#include <float.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define DBOR 1

typedef struct view_t
{
  struct camera_t cam[2]; // left and right eye cameras
  int num_cams;           // number of cameras in this view
  float eye_dist;         // do a stereo render
  float active_camid;     // camera currently active in gui (float for display params)

  uint64_t width;         // dimensions of frame buffer
  uint64_t height;        // and output renders

  uint64_t overlays;      // number of progressions
  float gain;             // extra gain for instance for mlt mean brightness

  double time_wallclock;  // render times
  double time_overlays;   // accumulated wallclock time per progression (no prepare etc)
  double time_user;

  double *stat_enery;     // path statistics
  uint64_t *stat_cnt;

  int num_fbs;            // number of frame buffers. can be != num_cams for gpt
  framebuffer_t *fb;      // colour frame buffers
  
  // dbor buffer cascade, currently not working for stereo
  int num_dbors;          // number of dbor buffers in cascade
  framebuffer_t *dbor;    // density based outlier rejection buffer cascade
  framebuffer_t *gcov;    // gaussian splat covariances per pixel per dbor buffer

  int lf_tile_size;       // tile size for light field camera
  float lf_scale;         // light field camera scale for directions (converting mm offsets to fraction of tile)

  int welch;                    // also record stuff needed for welch tests?
  double fb_welch_sample_count; // basically samples per pixel/3
  double *fb_welch_squared;     // records sum of (square of sum of 32x32 pixel block per frame) - 1 sample in welch is the sum of a 32x32 pixel block with 3 samples per pixel
  double *fb_welch_sum;         // records sum of welch samples. This is different from the frame buffer as it updates only every 3 samples per pixel
  double *fb_welch_tmp;         // records current 3 frames. Each 3 frames this is added squared to fb_welch_squared, and added to fb_welch_sum

  view_move_t moving;
}
view_t;

static const float view_full_frame_width = 0.35f; // [mm]
static const float view_f_stop[] = {
  0.5, 0.7, 1.0, 1.4, 2, 2.8, 4, 5.6, 8, 11, 16, 22, 32, 45, 64, 90, 128
};  // [f/#]
static const int view_num_f_stop = sizeof(view_f_stop)/sizeof(view_f_stop[0]);
static const float view_exposure_time[] =
{
  60.0f, 30.0f, 15.0f, 8.0f, 4.0f, 2.0f, 1.0f, 0.5f, 1.0/4.0f, 1.0/8.0f,
  1.0/15.0f, 1.0/30.0f, 1.0/60.0f, 1.0/125.0f, 1.0/250.0f, 1.0/500.0f, 1.0/1000.0f, 1.0/2000.0f, 1.0/4000.0f, 1.0/8000.0f
}; // [s]
static const int view_num_exposure_time = sizeof(view_exposure_time)/sizeof(view_exposure_time[0]);

float view_tv2time(int tv)
{
  assert(tv >= 0 && tv < view_num_exposure_time);
  return view_exposure_time[tv];
}

float view_av2fstop(int av)
{
  assert(av >= 0 && av < view_num_f_stop);
  return view_f_stop[av];
}

uint64_t view_width()
{
  return rt.view->width;
}

uint64_t view_height()
{
  return rt.view->height;
}

uint64_t view_overlays()
{
  return rt.view->overlays;
}

void view_clear_frame(view_t *r)
{
  for(int c=0;c<r->num_fbs;c++)
    fb_clear(r->fb+c);
  if(r->num_dbors > 1) for(int c=0;c<r->num_dbors;c++)
  {
    fb_clear(r->dbor+c);
    fb_clear(r->gcov+c);
  }

  r->gain = 1.0f;
  r->overlays = 0;
  r->time_wallclock = common_time_wallclock();
  r->time_overlays = 0.0;
  r->time_user = common_time_user();
  memset(r->stat_enery, 0, sizeof(double)  *(PATHSPACE_MAX_VERTS+1)*rt.num_threads);
  memset(r->stat_cnt,   0, sizeof(uint64_t)*(PATHSPACE_MAX_VERTS+1)*rt.num_threads);
}

void view_clear()
{
  view_clear_frame(rt.view);
  sampler_clear(rt.sampler);
  pointsampler_clear();
  render_clear();
}

static inline int parse_Av(const char *fstop)
{
  float fnum;
  sscanf(fstop, "f/%f", &fnum);
  for(int i=0;i<view_num_f_stop-1;i++)
    if(view_f_stop[i+1] > fnum) return i;
  return 2;
}

static void cam_init(const view_t *v, camera_t *c)
{
  strncpy(c->magic, "CCAM", 4);
  c->version = 1;

  // default cam focuses and looks at about the origin (z-up)
  // similar to default blender viewport.
  c->focus = 10.f;
  float axis[3] = {0, 1.0f/sqrtf(2.0f), 1.0f/sqrtf(2.0f)};
  quaternion_init(&(c->orient), M_PI, axis);
  c->orient_t1 = c->orient;
  c->pos_t1[0] = c->pos[0] = 0.0f;
  c->pos_t1[1] = c->pos[1] = -10.0f;
  c->pos_t1[2] = c->pos[2] = 1.0f;
  c->speed = 1.0f;

  c->aperture_value = 9;
  c->exposure_value = 13;
  c->crop_factor = 1.0f;
  c->focal_length = 0.5f;
  c->iso = 100.0f;

  if(v->width > v->height)
  {
    c->film_width = view_full_frame_width / c->crop_factor;
    c->film_height = view_height()/(float)view_width() * c->film_width;
  }
  else
  { // we're holding the camera in portrait mode:
    c->film_height = view_full_frame_width / c->crop_factor;
    c->film_width  = view_width()/(float)view_height() * c->film_height;
  }
}

void view_set_exposure_gain(const float gain)
{
  rt.view->gain = gain;
}

// helper to init the right eye from the left + eye distance
static void cam_init_stereo(
    const camera_t * const left,
    camera_t* right,
    const float eye_dist)
{
  if(eye_dist <= 0.0) return;
  // copy left eye values
  *right = *left;

  // check for motion blur
  char mb = (
      left->orient.x[0] == left->orient_t1.x[0] &&
      left->orient.x[1] == left->orient_t1.x[1] &&
      left->orient.x[2] == left->orient_t1.x[2] &&
      left->orient.w    == left->orient_t1.w &&
      left->pos[0] == left->pos_t1[0] &&
      left->pos[1] == left->pos_t1[1] &&
      left->pos[2] == left->pos_t1[2]
      ) ? 0 : 1;

  // we want to rotate the camera from the focal distance so that it
  // has the same focus and a fixed distance from one another. We just
  // find the angle between them and rotate
  const float alpha = 2.0f*asinf(0.5f * eye_dist/left->focus);

  for(int m=0;m<=mb;m++)
  {
    const quaternion_t *lefto = m ? &left->orient_t1 : &left->orient;
    quaternion_t *righto = m ? &right->orient_t1 : &right->orient;
    const float *leftp = m ? left->pos_t1 : left->pos;
    float *rightp = m ? right->pos_t1 : right->pos;

    float up[3] = {0, 1, 0};
    quaternion_transform(lefto, up);
    normalise(up);

    float fw[3] = {0, 0, 1};
    quaternion_transform(lefto, fw);
    normalise(fw);

    float focus_p[3];
    for(int i=0; i<3; i++) focus_p[i] = leftp[i] + fw[i]*left->focus;

    float rot_p[3];
    for(int i=0; i<3; i++) rot_p[i] = leftp[i] - focus_p[i];

    quaternion_t rot;
    quaternion_init(&rot, alpha, up);
    quaternion_transform(&rot, rot_p);

    for(int i=0; i<3; i++) rightp[i] = rot_p[i] + focus_p[i];

    quaternion_mult_fleft(&rot, righto);

    if(!mb)
    {
      right->orient_t1 = right->orient;
      for(int i=0; i<3; i++) right->pos_t1[i] = right->pos[i];
    }
  }
}

view_t *view_init()
{
  camera_global_init();
  view_t *v = (view_t *)malloc(sizeof(view_t));
  memset(v, 0, sizeof(*v));
  rt.view = v; // cheat a bit with initialisation
  v->num_cams = 1;
  v->num_fbs = 1;
  v->num_dbors = 0; // default no dbor
  v->eye_dist = 0.0f; // default off
  v->active_camid = 0;
  v->moving = 0;

  // default to 1k 16:9 render
  v->width = 1024;
  v->height = 576;

  // parse command line args
  const char *cam_file = 0;      // camera file to load
  float ISO = -1.0f;             // iso speed override
  int32_t Av = -3;               // aperture value Av override
  float focusdist = -1.0f;       // focus distance override
  int32_t Tv = -1;               // exposure time value Tv override
  float refocus_x = -1, refocus_y = -1;
  int fb_retain = 0;
  v->lf_tile_size = 0;           // light field stuff: pixel tile size
  v->lf_scale = 0.1f;            // scale converting directions to relative position in tile
  
  for(int i=0;i<rt.argc;i++)
  {
    if     ((strcmp(rt.argv[i], "-c") == 0) && (rt.argc > i+1)) cam_file = rt.argv[++i];
    else if((strcmp(rt.argv[i], "-w") == 0) && (rt.argc > i+1)) v->width = atol(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "-h") == 0) && (rt.argc > i+1)) v->height = atol(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--Av")    == 0) && (rt.argc > i+1)) Av = parse_Av(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--focus") == 0) && (rt.argc > i+1)) focusdist = atof(rt.argv[++i])*10.0f;
    else if((strcmp(rt.argv[i], "--Tv")    == 0) && (rt.argc > i+1)) Tv = atol(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--iso")   == 0) && (rt.argc > i+1)) ISO = atof(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--refocus") == 0) && (rt.argc > i+2)) { refocus_x = atof(rt.argv[++i]); refocus_y = atof(rt.argv[++i]); }
    else if((strcmp(rt.argv[i], "--eye-dist") == 0) && (rt.argc > i+1)) { v->eye_dist = atof(rt.argv[++i]); v->num_cams = 2; }
    else if( strcmp(rt.argv[i], "--gpt") == 0) { v->num_fbs = v->num_cams*3; }
    else if( strcmp(rt.argv[i], "--retain-framebuffer") == 0) { fb_retain = 1; }
    // XXX else if( strcmp(rt.argv[i], "--welch") == 0) { v->welch = 1; }
    else if((strcmp(rt.argv[i], "--lf-tile-size") == 0) && (rt.argc > i+1)) v->lf_tile_size = atol(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--lf-scale") == 0) && (rt.argc > i+1)) v->lf_scale = atof(rt.argv[++i]);
    else if((strcmp(rt.argv[i], "--dbor") == 0) && (rt.argc > i+1)) { v->num_dbors = atol(rt.argv[++i]); v->num_dbors = CLAMP(v->num_dbors, 0, 20); }
  }

  // fullfill 32-alignment for tiles.
  while(v->width  & 0x1f) v->width++;
  while(v->height & 0x1f) v->height++;

  // init camera format, if given:
  if(cam_file && strstr(cam_file, "%"))
    strncpy(rt.camformat, cam_file, 256);
  else snprintf(rt.camformat, 256, "%s%%02d.cam", rt.basename);

  // init default camera
  for(int k=0;k<v->num_cams;k++)
    cam_init(v, v->cam);

  // load camera, if any (will update display and clear view, using the thread
  // pool. depends on these to be inited)
  if(!cam_file || view_cam_read(cam_file))
  { // try basename01.cam
    char filename[512] = {0};
    snprintf(filename, sizeof(filename), rt.camformat, 1);
    view_cam_read(filename);
  }

  // init welch test variance buffer
  v->fb_welch_squared = 0;
  v->fb_welch_sum = 0;
  v->fb_welch_tmp = 0;
  if(v->welch)
  {
    int w_wd = v->width / 32;
    int w_ht = v->height / 32;
    v->fb_welch_squared = malloc(sizeof(double)*3*w_wd*w_ht);
    v->fb_welch_sum = malloc(sizeof(double)*3*w_wd*w_ht);
    v->fb_welch_tmp = malloc(sizeof(double)*3*w_wd*w_ht);
  }

  // init frame buffers
  char filename[1024];
  v->fb = malloc(sizeof(framebuffer_t)*v->num_fbs);
  for(int c=0;c<v->num_fbs;c++)
  {
    snprintf(filename, sizeof(filename), "%s_%s_fb%02d.fb", rt.basename, rt.output_filename, c);
    fb_init(v->fb + c, v->width, v->height, 3, filename);
    v->fb[c].retain = fb_retain;
  }
  v->dbor = 0;
  if(v->num_dbors > 1)
  {
    v->dbor = malloc(sizeof(framebuffer_t)*v->num_dbors);
    v->gcov = malloc(sizeof(framebuffer_t)*v->num_dbors);
    for(int c=0;c<v->num_dbors;c++)
    {
      snprintf(filename, sizeof(filename), "%s_%s_dbor%02d.fb", rt.basename, rt.output_filename, c);
      fb_init(v->dbor + c, v->width, v->height, 3, filename);
      v->dbor[c].retain = fb_retain;
      fprintf(stderr, "[view] dbor cascade buffer [%02d] %s\n", c, filename);
      snprintf(filename, sizeof(filename), "%s_%s_gcov%02d.fb", rt.basename, rt.output_filename, c);
      fb_init(v->gcov + c, v->width, v->height, 3, filename);
      v->gcov[c].retain = fb_retain;
    }
  }

  v->stat_enery = (double   *)malloc((PATHSPACE_MAX_VERTS+1)*rt.num_threads*sizeof(double));
  v->stat_cnt   = (uint64_t *)malloc((PATHSPACE_MAX_VERTS+1)*rt.num_threads*sizeof(uint64_t));
  view_clear_frame(v);


  camera_t *cam = v->cam + 0;
  if(focusdist != -1.0)
  {
    camera_set_focus(cam, focusdist);
    printf("[view] focus distance %fm\n", .1f*focusdist);
  }
  if(ISO > 0)
  {
    cam->iso = ISO;
    printf("[view] iso %f\n", ISO);
  }
  if(Tv >= 0 && Tv <= view_num_exposure_time)
  {
    cam->exposure_value = Tv;
    if(cam->exposure_value > 6)
      printf("[view] Tv : 1/%.0f\n", roundf(1.0/view_exposure_time[cam->exposure_value]));
    else
      printf("[view] Tv : %.1f''\n", view_exposure_time[cam->exposure_value]);
  }
  if(Av >= 0 && Av <= view_num_f_stop)
  {
    cam->aperture_value = Av;
    printf("[view] Av : f/%.1f\n", view_f_stop[cam->aperture_value]);
  }
  if(refocus_x >= 0 && refocus_y >= 0)
  {
    view_set_focus(refocus_x, refocus_y, 1);
  }
  cam_init_stereo(cam, cam+1, v->eye_dist);

  rt.display = display_open("corona-13", view_width(), view_height());
  display_control_add(rt.display, "[cam] id ", &rt.view->active_camid, 0, rt.view->num_cams-1, 1, 0, 0);

  return v;
}

void view_cleanup(view_t *v)
{
  free(v->stat_enery);
  free(v->stat_cnt);
  for(int c=0;c<v->num_fbs;c++)
    fb_cleanup(v->fb+c);
  if(v->num_dbors > 1) for(int c=0;c<v->num_dbors;c++)
  {
    fb_cleanup(v->dbor+c);
    fb_cleanup(v->gcov+c);
  }
  if(v->welch)
  {
    free(v->fb_welch_squared);
    free(v->fb_welch_sum);
    free(v->fb_welch_tmp);
  }
  free(v->fb);
  free(v->dbor);
  free(v->gcov);
  free(v);
}

const camera_t *view_get_camera(const path_t *path)
{
  const int camid = path->sensor.camid;
  assert(camid >= 0);
  assert(camid < rt.view->num_cams);
  return rt.view->cam + camid;
}

int view_sample_camid(const float rand)
{
  return CLAMP(rt.view->num_cams * rand, 0, rt.view->num_cams-1);
}

float view_pdf_camid(const int camid)
{
  if(camid < 0 || camid >= rt.view->num_cams) return 0.0f;
  return 1.0f/rt.view->num_cams;
}

void view_splat_gaussian(const path_t *path, const float mu[2], const float S[4], const float col[3])
{
  // const int tid = common_get_threadid();
  // rt.view->stat_enery[(PATHSPACE_MAX_VERTS+1)*tid + path->length] += value;
  // rt.view->stat_cnt  [(PATHSPACE_MAX_VERTS+1)*tid + path->length] ++;

  // TODO: call deferred splat as the others are doing?
  // TODO: splat to dbor and gcov buffers!
  // TODO: how to normalise these properly?

  const int cid = path->sensor.camid;
  filter_gaussian_splat(rt.view->fb+cid,
    path->sensor.pixel_i, path->sensor.pixel_j, mu, S, col);

  // TODO: bring dbor back here and let the splatting infrastructure take the framebuffer, too!
}

void view_splat(const path_t *path, const mf_t value)
{
  if(!mf_any(mf_gt(value, mf_set1(0.0f)))) return; // fast path for == 0.0
  if(!mf_all(mf_lt(value, mf_set1(FLT_MAX)))) return; // reject inf
  if(!mf_all(mf_eq(value, value))) return; // filter NaN
  float col[3];
  view_deferred_splat(path, value, col);
  view_splat_col(path, col);
}

void view_deferred_splat(const path_t *path, const mf_t value, float *col)
{
  col[0] = col[1] = col[2] = 0.0f;
  if(!mf_any(mf_gt(value, mf_set1(0.0f)))) return;
  const int tid = common_get_threadid();
  rt.view->stat_enery[(PATHSPACE_MAX_VERTS+1)*tid + path->length] += mf_hsum(value);
  rt.view->stat_cnt  [(PATHSPACE_MAX_VERTS+1)*tid + path->length] ++;
  spectrum_p_to_camera(path->lambda, value, col);
}

void view_splat_col(const path_t *path, const float *col)
{
  const int cid = path->sensor.camid;
  if(rt.view->lf_tile_size)
  { // light field camera, don't filter:
    const int ts = rt.view->lf_tile_size;
    const float psf_scale = rt.view->lf_scale;
    const float rad = ts/(2.0 * psf_scale);
    int i = path->sensor.pixel_i / ts, j = path->sensor.pixel_j / ts;
    i *= ts; j *= ts;
    i += ts/2 + path->sensor.pixel_dir_i * rad;
    j += ts/2 + path->sensor.pixel_dir_j * rad;

    if(i >= 0 && i < view_width() && j >= 0 && j < view_height())
      filter_box_splat(rt.view->fb+cid, i, j, col);
  }
  else
  { // vanilla splat
    filter_splat(rt.view->fb+cid,
        path->sensor.pixel_i, path->sensor.pixel_j, col);
  }

  if(rt.view->num_dbors > 1)
  {
    const float value = (col[0] + col[1] + col[2])/3.f;
    if (value > 0)
    {
      float coll[3], colu[3];
      const float logval = MAX(0, log2f(value));
      const int l = CLAMP((int)logval, 0, rt.view->num_dbors-1);
      const int u = l+1;
      const float lv = CLAMP(value < 1.f ? 1.f : (((1<<l)/value)-0.5f)/0.5f, 0.0f, 1.0f);
      const float uv = 1.f - lv;
      if(l == rt.view->num_dbors-1)
        for (int k = 0; k < 3; ++k)
          coll[k] = col[k];
      else for (int k = 0; k < 3; ++k)
      {
        coll[k] = lv * col[k];
        colu[k] = uv * col[k];
      }
      filter_blackmanharris_splat(rt.view->dbor+l,
          path->sensor.pixel_i, path->sensor.pixel_j, coll);
      if (u < rt.view->num_dbors)
        filter_blackmanharris_splat(rt.view->dbor+u,
            path->sensor.pixel_i, path->sensor.pixel_j, colu);
    }
  }

#if 0
  // TODO: bring this feature back in alisa's version (and doubles)
  if(rt.view->welch)
  {
    // fill tmp welch buffer with sum over 32x32 pixels and 3 samples per pixel
    const int w_wd = rt.view->width / 32;
    const int w_ht = rt.view->height / 32;
    const int i = path->sensor.pixel_i / 32.0f;
    const int j = path->sensor.pixel_j / 32.0f;
    if(i >= 0 && i < w_wd && j >= 0 && j < w_ht)
    {
      double *p = rt.view->fb_welch_tmp + 3*(i+w_wd*j);
      for(int k=0;k<3;k++) common_atomic_add64(p+k, (double)(col[k]));
    }
  }
#endif
}

// output render to file:
void view_write_images(const char *suffix)
{
  char filename[1024];
  pointsampler_finalize(rt.pointsampler);
  for(int fid=0;fid<rt.view->num_fbs;fid++)
  {
    snprintf(filename, sizeof(filename), "%s%s_fb%02d.pfm", rt.basename, suffix, fid);
    fb_export(rt.view->fb+fid, filename, 0, 3);
    common_write_sidecar(filename);
  }
  if(rt.view->num_dbors > 1) for(int fid=0;fid<rt.view->num_dbors;fid++)
  {
    snprintf(filename, sizeof(filename), "%s%s_dbor%02d.pfm", rt.basename, suffix, fid);
    fb_export(rt.view->dbor+fid, filename, 0, 3);
    snprintf(filename, sizeof(filename), "%s%s_gcov%02d.pfm", rt.basename, suffix, fid);
    fb_export(rt.view->gcov+fid, filename, 0, 3);
  }
#if 0
  if(rt.view->welch)
  {
    printf("welch square0: %f, welch sum0: %f\n",
        rt.view->fb_welch_squared[0],
        rt.view->fb_welch_sum[0]);

    // write sum of squares
    snprintf(filename, sizeof(filename), "%s%s_welch2", rt.basename, suffix);
    const camera_t *cam = rt.view->cam + 0;
    const float inv_o = rt.view->gain * cam->iso/100.0f;
    FILE *f1 = fopen(filename, "wb");
    if(f1)
    {
      const uint64_t w_wd = rt.view->width / 32, w_ht = rt.view->height/32;
      for(uint64_t k = 0; k < w_wd*w_ht; k++)
      {
        double w[3] = {
          rt.view->fb_welch_squared[3*k],
          rt.view->fb_welch_squared[3*k+1],
          rt.view->fb_welch_squared[3*k+2]
        };
        for(int i=0;i<3;i++) w[i] *= inv_o*inv_o; // scale sum of squares [0-2]
        fwrite(w, sizeof(double), 3, f1);
      }
      double N[1] = {rt.view->fb_welch_sample_count}; // write number of welch samples
      fwrite(N, sizeof(double), 1, f1);
      double ovl[1] = {rt.view->overlays}; // write number of samples per pixel
      fwrite(ovl, sizeof(double), 1, f1);
      fclose(f1);
    }

    // write sum
    snprintf(filename, sizeof(filename), "%s%s_welch1", rt.basename, suffix);
    FILE *f2 = fopen(filename, "wb");
    if(f2)
    {
      const uint64_t w_wd = rt.view->width / 32, w_ht = rt.view->height/32;
      for(uint64_t k = 0; k < w_wd*w_ht; k++)
      {
        double w[3] = {
          rt.view->fb_welch_sum[3*k],
          rt.view->fb_welch_sum[3*k+1],
          rt.view->fb_welch_sum[3*k+2]
        };
        for(int i=0;i<3;i++) w[i] *= inv_o; // scale sum
        fwrite(w, sizeof(double), 3, f2);
      }
      double N[1] = {rt.view->fb_welch_sample_count};
      fwrite(N, sizeof(double), 1, f2);
      double ovl[1] = {rt.view->overlays};
      fwrite(ovl, sizeof(double), 1, f2);
      fclose(f2);
    }
  }
#endif
}

static void *work_sample(void *arg)
{
  threads_t *job = (threads_t *)arg;
  while(1)
  {
    uint64_t i =  __sync_fetch_and_add(&job->counter, 1);
    if(i >= job->end) return 0;
    // fprintf(stderr, "[worker %03d] processing sample %lu\n", rt_tls.tid, i);
    render_sample_path(i);
  }
}

void view_render()
{
  threads_t *t = rt.threads;

  const double time_begin = common_time_wallclock();
  // step sample counter globally in 1spp intervals:
  const uint64_t end = t->counter + (uint64_t)rt.batch_frames * (uint64_t)(rt.view->width * rt.view->height);
  t->counter = t->end;
  t->end = end;
  // do this only now so they know how many samples the threads use
  lights_prepare_frame();
  sampler_prepare_frame(rt.sampler);
  pointsampler_prepare_frame(rt.pointsampler);
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(t->task + k, &t->pool, work_sample, t);
  pthread_pool_wait(&t->pool);

  // accumulate last paths in markov chain and update gain
  pointsampler_finalize(rt.pointsampler);

  // update progression count in frame buffer files:
  rt.view->overlays += rt.batch_frames;
  for(int c=0;c<rt.view->num_fbs;c++)
  {
    const int cid = c / (rt.view->num_fbs/rt.view->num_cams);
    const camera_t *cam = rt.view->cam + cid;
    rt.view->fb[c].header->gain = rt.view->gain * cam->iso / (100.0f * rt.view->overlays);
  }

  for(int c=0;c<rt.view->num_dbors;c++)
  {
    const int cid = 0;
    const camera_t *cam = rt.view->cam + cid;
    rt.view->dbor[c].header->gain = rt.view->gain * cam->iso / (100.0f * rt.view->overlays);
    rt.view->gcov[c].header->gain = rt.view->gain * cam->iso / (100.0f * rt.view->overlays);
  }

  if (rt.view->welch)
  {
    if (rt.view->overlays % 3 == 0)
    {
      // every 3 frames: accumulate welch buffers
      const int w_wd = rt.view->width / 32;
      const int w_ht = rt.view->height / 32;
      printf("clear welch tmp %f %f %f\n", rt.view->fb_welch_tmp[0], rt.view->fb_welch_squared[0], rt.view->fb_welch_sample_count);
      for (int i=0;i<w_wd*w_ht*3;i++)
      {
        // add square
        common_atomic_add64(rt.view->fb_welch_squared+i, rt.view->fb_welch_tmp[i] * rt.view->fb_welch_tmp[i]);
        // add value
        common_atomic_add64(rt.view->fb_welch_sum+i, rt.view->fb_welch_tmp[i]);
        rt.view->fb_welch_tmp[i] = 0.0;
      }
      rt.view->fb_welch_sample_count += 1.0;
    }
  }

  const double time_end = common_time_wallclock();
  rt.view->time_overlays += time_end - time_begin;

  // trigger display update
  display_update(
      rt.display,
      rt.view->fb[view_get_active_camid()].fb,
      rt.view->fb[view_get_active_camid()].header->gain);
}

// gui interaction, movement:
void view_set_active_camid(int camid)
{
  assert(camid >= 0);
  assert(camid < rt.view->num_cams);
  rt.view->active_camid = camid;
}

int view_get_active_camid()
{
  return rt.view->active_camid;
}

static inline void _view_print_time(FILE *fd, double time)
{
  int seconds = (int)time;
  int milli = (time - seconds)*1000;
  int mins = seconds / 60.0;
  seconds -= mins * 60;
  int hours = mins / 60.0;
  mins -= hours * 60;
  int days = hours / 24.0;
  hours -= days * 24;
  if(days > 0) fprintf(fd, "%d days ", days);
  if(days > 0 || hours > 0) fprintf(fd, "%dh ", hours);
  if(days > 0 || hours > 0 || mins > 0) fprintf(fd, "%02d:", mins);
  fprintf(fd, "%02d.%02d", seconds, milli/10);
}

void view_print_info(FILE *fd)
{
  const double wallclock = common_time_wallclock() - rt.view->time_wallclock;
  const double user = common_time_user() - rt.view->time_user;
  fprintf(fd, "view     : samples per pixel: %lu (%.2f s/prog) max path vertices %d\n", rt.view->overlays, rt.view->time_overlays/rt.view->overlays, PATHSPACE_MAX_VERTS);
  fprintf(fd, "           res %lux%lu\n", rt.view->width, rt.view->height);
  fprintf(fd, "           elapsed wallclock prog %.2fs (", rt.view->time_overlays);
  _view_print_time(fd, rt.view->time_overlays);
  fprintf(fd, "), total %.2fs (", wallclock);
  _view_print_time(fd, wallclock);
  fprintf(fd, "), user %.2fs (", user);
  _view_print_time(fd, user);
  fprintf(fd, ")\n");

  fprintf(fd, "           active cam %d\n", (int)rt.view->active_camid);
  int fbs_p_cam = rt.view->num_fbs/rt.view->num_cams;
  for(int c=0;c<rt.view->num_cams;c++)
  {
    camera_print_info(rt.view->cam + c, fd);
    const float inv_o = rt.view->gain * rt.view->cam[c].iso/100.0f * 1.0f/rt.view->overlays;
    const float inv_px = inv_o/(rt.view->width*rt.view->height);
    float sumr = 0, sumg = 0, sumb = 0;
    // #pragma omp parallel for default(none) shared(rt) reduction(+:sumr) reduction(+:sumg) reduction(+:sumb) schedule(static)
    for(int k=0;k<3*rt.view->width*rt.view->height;)
    {
      sumr += inv_px*rt.view->fb[c*fbs_p_cam].fb[k++];
      sumg += inv_px*rt.view->fb[c*fbs_p_cam].fb[k++];
      sumb += inv_px*rt.view->fb[c*fbs_p_cam].fb[k++];
    }
    fprintf(fd, "           cam %d average image intensity (rgb): (%f %f %f)\n", c, sumr, sumg, sumb);
  }

  const int lines = 3;
  double max = 0.0;
  double base = 64.0;
  double p[PATHSPACE_MAX_VERTS+1] = {0.0};
  for(int k=2;k<=PATHSPACE_MAX_VERTS;k++) for(int t=0;t<rt.num_threads;t++)
  {
    if(rt.view->stat_cnt[(PATHSPACE_MAX_VERTS+1)*t + k])
      p[k] += rt.view->stat_enery[(PATHSPACE_MAX_VERTS+1)*t + k];
    max = MAX(max, p[k]);
  }
  fprintf(fd, "           ");
  for(int k=2;k<=PATHSPACE_MAX_VERTS;k++) if(k % 10 == 0) fprintf(fd, "%d", (k/10)%10); else if(k%5==0) fprintf(fd, "|"); else fprintf(fd, ".");
  fprintf(fd, "\n");
  for(int h=lines-1;h>=0;h--)
  {
    fprintf(fd, "           ");
    for(int k=2;k<=PATHSPACE_MAX_VERTS;k++)
    {
      float level = p[k] / max;
      level = log(1.0 + (base-1.0)*level)/log(base);
      float fill = level*lines - h;
      if(fill <= 0) fprintf(fd, " ");
      else if(fill <= 1./8.) fprintf(fd, "\u2581");
      else if(fill <= 2./8.) fprintf(fd, "\u2582");
      else if(fill <= 3./8.) fprintf(fd, "\u2583");
      else if(fill <= 4./8.) fprintf(fd, "\u2584");
      else if(fill <= 5./8.) fprintf(fd, "\u2585");
      else if(fill <= 6./8.) fprintf(fd, "\u2586");
      else if(fill <= 7./8.) fprintf(fd, "\u2587");
      else /*if(fill <= 8./8.)*/ fprintf(fd, "\u2588");
    }
    fprintf(fd, "\n");
  }

  render_print_info(fd);
  filter_print_info(fd);
}

void view_set_focus(float x, float y, int verbose)
{
  path_t path;
  path_init(&path, lrand48(), 0);
  path.sensor.pixel_i = x;
  path.sensor.pixel_j = y;
  path.sensor.pixel_set = 1;
  path.lambda = mf_set1(550.0f);
  if(path_extend(&path)) return;
  if(path.e[1].dist < FLT_MAX)
  {
    camera_t *cam = rt.view->cam + 0;
    camera_set_focus(cam, path.e[1].dist);
    cam_init_stereo(cam, cam+1, rt.view->eye_dist);
    if(verbose)
    {
      printf("camera focus        : at pixel (%d, %d) dist %f m\n", (int)x, (int)y, path.e[1].dist*0.1f);
      if(!primid_invalid(path.v[1].hit.prim))
          printf("shape under mouse   : %s\n", rt.prims->shape[path.v[1].hit.prim.shapeid].name);
    }
    view_clear();
  }
}

void view_print_usage()
{
  fprintf(stderr, "  [-w width -h height]                    image dimensions\n"
                  "  [--focus focusdistance(m)]              focus distance in metres\n"
                  "  [--Av Av]                               aperture e.g. f/4.0\n" 
                  "  [--Tv Tv]                               exposure time value\n"
                  "  [--iso ISO]                             film back iso\n"
                  "  [--refocus <x> <y>]                     after loading cam, focus on this pixel\n"
                  "  [--eye-dist <x>]                        adjust eye distance, 0 for mono render\n"
                  "  [--retain-framebuffer]                  don't unlink frame buffer on shutdown\n"
                  "  [--lf-tile-size <x>]                    light field tile size in pixels\n"
                  "  [--lf-scale <x>]                        scale direction of light field to relative tile position\n"
                  "                                          (this only depends on the lens, not on resolution)\n"
                  "  [--dbor <x>]                            use x levels of density based outlier rejection buffers\n"
                  "  [--welch]                               output buffer for welch test\n");
}


void view_cam_display_info()
{
  camera_display_info(rt.view->cam+(int)rt.view->active_camid);
}

mf_t view_cam_sample(path_t *p)
{
  // sample camera, and compensate for uniform pdf = 1/num_cams
  if(!p->sensor.camid_set)
    p->sensor.camid = view_sample_camid(pointsampler(p, s_dim_camid));
  return mf_mul(mf_set1(rt.view->num_cams), camera_sample(rt.view->cam+p->sensor.camid, p));
}

mf_t view_cam_pdf(const path_t *p, int v)
{
  return camera_pdf(rt.view->cam+p->sensor.camid, p, v);
}

mf_t view_cam_connect(path_t *p)
{
  return camera_connect(rt.view->cam+p->sensor.camid, p);
}

mf_t view_cam_pdf_connect(const path_t *p, int v)
{
  return camera_pdf_connect(rt.view->cam+p->sensor.camid, p, v);
}

mf_t view_cam_mutate_aperture(path_t *p, const float r1, const float r2, const float step)
{
  return camera_mutate_aperture(rt.view->cam+p->sensor.camid, p, r1, r2, step);
}

mf_t view_cam_pdf_mutate_aperture(const path_t *curr, const path_t *tent, const float step)
{
  return camera_pdf_mutate_aperture(rt.view->cam+curr->sensor.camid, curr, tent, step);
}

mf_t view_cam_eval(const path_t *p)
{
  return camera_eval(rt.view->cam+p->sensor.camid, p);
}

float view_sample_time(const path_t *p, const float rand)
{
  const camera_t *c = view_get_camera(p);
  // time interval normalised to 1/30 ~ 24 fps @270deg shutter (old style cinema)
  // we don't support exposure longer than that.
  // 1/30 can be exactly represented as shutter time as exposure_value.
  const float maxexp = 1.0f/30.0f; // max exposure time in seconds
  // exposure time in seconds:
  const float exposure_time = view_tv2time(c->exposure_value);
  return rand * fminf(1.0f, exposure_time/maxexp);
}

float view_pdf_time(const path_t *p, const float time)
{
  const camera_t *c = view_get_camera(p);
  const float maxexp = 1.0f/30.0f;
  const float exposure_time = view_tv2time(c->exposure_value);
  // irrelevant as long as its constant per image:
  return 1./fminf(maxexp, exposure_time);
}

// init the tangent frame pointed to by *hit for the given time.
void view_cam_init_frame(const path_t *p, hit_t *hit)
{
  const camera_t *c = rt.view->cam + (p ? p->sensor.camid : 0);
  const float time = p ? p->time : 0.0f;
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

int view_cam_write(const char *const filename)
{
  camera_t *cam = rt.view->cam + (int)rt.view->active_camid;
  return camera_write(cam, filename);
}

int view_cam_write_readable(const char *const filename_readable)
{
  camera_t *cam = rt.view->cam + (int)rt.view->active_camid;
  return camera_write_humanly_readable(cam, filename_readable);
}

int view_cam_read(const char *const filename)
{
  camera_t *cam = rt.view->cam + 0;
  cam_init(rt.view, cam);
  if(camera_read(cam, filename)) return 1;
  if(cam->exposure_value < 0 || cam->exposure_value > view_num_exposure_time) cam->exposure_value = 13;
  if(rt.view->width > rt.view->height)
  {
    cam->film_width = view_full_frame_width / cam->crop_factor;
    cam->film_height = view_height()/(float)view_width() * cam->film_width;
  }
  else
  { // we're holding the camera in portrait mode:
    cam->film_height = view_full_frame_width / cam->crop_factor;
    cam->film_width  = view_width()/(float)view_height() * cam->film_height;
  }
  if(cam->iso < 1 || cam->iso > 409600) cam->iso = 100;
  camera_set_focus(cam, cam->focus);
  cam_init_stereo(cam, cam+1, rt.view->eye_dist);
  return 0;
}


// user interaction

void view_active_cam_rotate(const float *axis, const float angle)
{
  camera_t *cam = rt.view->cam + 0;
  quaternion_t rot;
  quaternion_init(&rot, angle, axis);
  quaternion_mult(&(cam->orient), &rot);
  quaternion_mult(&(cam->orient_t1), &rot);
  cam_init_stereo(cam, cam+1, rt.view->eye_dist);
  view_clear();
}

void view_active_cam_move()
{
  camera_t *c = rt.view->cam + 0;
  view_move_t m = rt.view->moving;
  if(!m) return;
  hit_t hit;
  view_cam_init_frame(0, &hit);
  if(m & s_view_move_fw) for(int k=0;k<3;k++) c->pos   [k] += c->speed*hit.n[k];
  if(m & s_view_move_bk) for(int k=0;k<3;k++) c->pos   [k] -= c->speed*hit.n[k];
  if(m & s_view_move_rg) for(int k=0;k<3;k++) c->pos   [k] -= c->speed*hit.a[k];
  if(m & s_view_move_lf) for(int k=0;k<3;k++) c->pos   [k] += c->speed*hit.a[k];
  if(m & s_view_move_up) for(int k=0;k<3;k++) c->pos   [k] += c->speed*hit.b[k];
  if(m & s_view_move_dn) for(int k=0;k<3;k++) c->pos   [k] -= c->speed*hit.b[k];
  if(m & s_view_move_fw) for(int k=0;k<3;k++) c->pos_t1[k] += c->speed*hit.n[k];
  if(m & s_view_move_bk) for(int k=0;k<3;k++) c->pos_t1[k] -= c->speed*hit.n[k];
  if(m & s_view_move_rg) for(int k=0;k<3;k++) c->pos_t1[k] -= c->speed*hit.a[k];
  if(m & s_view_move_lf) for(int k=0;k<3;k++) c->pos_t1[k] += c->speed*hit.a[k];
  if(m & s_view_move_up) for(int k=0;k<3;k++) c->pos_t1[k] += c->speed*hit.b[k];
  if(m & s_view_move_dn) for(int k=0;k<3;k++) c->pos_t1[k] -= c->speed*hit.b[k];
  cam_init_stereo(c, c+1, rt.view->eye_dist);
  view_clear(); // moved, need to invalidate frame buffer
}

void view_move_begin(view_move_t m)
{
  rt.view->moving |= m;
}

void view_move_end(view_move_t m)
{
  rt.view->moving &= ~m;
}

void view_ctl(view_ctl_t c)
{
  camera_t *cam = rt.view->cam + 0;
  switch(c)
  {
    case s_view_ctl_av_dn:
      // adjust f-stop
      if(cam->aperture_value < view_num_f_stop) cam->aperture_value++;
      break;
    case s_view_ctl_av_up:
      if(cam->aperture_value > 0) cam->aperture_value--;
      break;
    case s_view_ctl_fl_dn:
      // adjust zoom
      cam->focal_length /= 6.0f/5.0f;
      break;
    case s_view_ctl_fl_up:
      cam->focal_length *= 6.0f/5.0f;
      break;
    case s_view_ctl_tv_dn:
      // adjust exposure time
      if(cam->exposure_value > 0) cam->exposure_value--;
      break;
    case s_view_ctl_tv_up:
      if(cam->exposure_value < view_num_exposure_time-1) cam->exposure_value++;
      break;
    case s_view_ctl_speed_dn:
      cam->speed /= 2.0f;
      printf("camera speed        : %g m/frame\n", cam->speed*0.1);
      break;
    case s_view_ctl_speed_up:
      cam->speed *= 2.0f;
      printf("camera speed        : %g m/frame\n", cam->speed*0.1);
      break;
    case s_view_ctl_iso_up:
      cam->iso *= 2.0f;
      break;
    case s_view_ctl_iso_dn:
      cam->iso /= 2.0f;
      break;
    default:
      break;
  }
  cam_init_stereo(cam, cam+1, rt.view->eye_dist);
  view_cam_display_info();
  if(c != s_view_ctl_iso_up && c != s_view_ctl_iso_dn)
    view_clear();
}
