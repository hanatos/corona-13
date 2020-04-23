/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "prims.h"
#include "accel.h"
#include "view.h"
#include "shader.h"
#include "threads.h"
#include "display.h"
#include "points.h"
#include "lights.h"
#include "camera.h"
#include "sampler.h"
#include "spectrum.h"
#include "rgb2spec.h"
#include "version.h"

#include <stdlib.h>
#include <unistd.h>

/* global setting */
rt_t rt;
__thread rt_tls_t rt_tls;

/* gui stuff */
int gui_display = 1;
#ifndef GUI_QWERTZ
#ifndef GUI_NEO2
#define gui_key_lft KeyA
#define gui_key_bck KeyO
#define gui_key_rgt KeyE
#define gui_key_fwd KeyComma
#define gui_key_dwn KeyJ
#define gui_key_sp  KeyPeriod
#define gui_key_sm  KeySemiColon
#else
#define gui_key_lft KeyU
#define gui_key_bck KeyI
#define gui_key_rgt KeyA
#define gui_key_fwd KeyV
#define gui_key_dwn KeyJ
#define gui_key_sp  KeyW
#define gui_key_sm  KeyO
#endif
#else
#define gui_key_lft KeyA
#define gui_key_bck KeyS
#define gui_key_rgt KeyD
#define gui_key_fwd KeyW
#define gui_key_dwn KeyX
#define gui_key_sp  KeyE
#define gui_key_sm  KeyY
#endif
#define gui_key_tp  KeyOne
#define gui_key_tm  KeyTwo
#define gui_key_ap  KeyThree
#define gui_key_am  KeyFour
#define gui_key_fp  KeyFive
#define gui_key_fm  KeySix
#define gui_key_ip  KeySeven
#define gui_key_im  KeyEight

/* event handlers */
void onKeyDown(keycode_t key)
{
  static int num = 1;
  switch (key)
  {
    case KeyZero: // clear frame buffer, but not pointsampler (guiding cache)
      view_clear_frame(rt.view);
      // pointsampler_stop_learning(rt.pointsampler);
      break;
    case KeyEscape:
      rt.quit = 1;
      break;
    case gui_key_fwd:
      view_move_begin(s_view_move_fw);
      break;
    case gui_key_bck:
      view_move_begin(s_view_move_bk);
      break;
    case gui_key_lft:
      view_move_begin(s_view_move_lf);
      break;
    case gui_key_rgt:
      view_move_begin(s_view_move_rg);
      break;
    case KeyC:
      {
        char filename[265];
        for(int k=0;k<10;k++)
        {
          sprintf(filename, "%s%02d.cam", rt.basename, num++);
          if(!view_cam_write(filename))
          {
            display_print(rt.display, 0, 0, "saved camera %d", num-1);
            break;
          }
        }
        break;
       }
    case KeyD:
      {
        char filename_readable[265];
        for(int k=0;k<10;k++)
        {
          sprintf(filename_readable, "%s%02d_cam.txt", rt.basename, num++);
          if(!view_cam_write_readable(filename_readable))
          {
            display_print(rt.display, 0, 0, "saved readable camera %d", num-1);
            break;
          }
        }
        break;
      }
    case KeyL:
      {
        char filename[265];
        //allow one retry (modulo) and one to jump over potentially missing 00.cam
        for(int i=0;i<3;i++)
        {
          sprintf(filename, rt.camformat, rt.anim_frame++);
          if(view_cam_read(filename) == 0)
          {
            display_print(rt.display, 0, 0, "loaded camera %d", rt.anim_frame-1);
            view_clear();
            break;
          } else if(rt.anim_frame>1) rt.anim_frame = 0;
        }
        break;
      }
    case KeyH:
      if(gui_display ^= 1) display_print(rt.display, 0, 0, "");
      else view_cam_display_info();
      break;
    case KeyP:
      rt.screenshot = 1;
      break;
    case KeySpace:
      view_move_begin(s_view_move_up);
      break;
    case gui_key_dwn:
      view_move_begin(s_view_move_dn);
      break;
    case gui_key_am:
      view_ctl(s_view_ctl_av_dn);
      break;
    case gui_key_ap:
      view_ctl(s_view_ctl_av_up);
      break;
    case gui_key_fp:
      view_ctl(s_view_ctl_fl_dn);
      break;
    case gui_key_fm:
      view_ctl(s_view_ctl_fl_up);
      break;
    case gui_key_tp:
      view_ctl(s_view_ctl_tv_dn);
      break;
    case gui_key_tm:
      view_ctl(s_view_ctl_tv_up);
      break;
    case gui_key_sm:
      view_ctl(s_view_ctl_speed_dn);
      break;
    case gui_key_sp:
      view_ctl(s_view_ctl_speed_up);
      break;
    case gui_key_ip:
      view_ctl(s_view_ctl_iso_up);
      break;
    case gui_key_im:
      view_ctl(s_view_ctl_iso_dn);
      break;
    default:
      break;
  }
}

void onKeyUp(keycode_t key)
{
  switch (key)
  {
    case gui_key_fwd:
      view_move_end(s_view_move_fw);
      break;
    case gui_key_bck:
      view_move_end(s_view_move_bk);
      break;
    case gui_key_lft:
      view_move_end(s_view_move_lf);
      break;
    case gui_key_rgt:
      view_move_end(s_view_move_rg);
      break;
    case KeySpace:
      view_move_end(s_view_move_up);
      break;
    case gui_key_dwn:
      view_move_end(s_view_move_dn);
      break;
    default:
      break;
  }
}

void onMouseButtonDown(mouse_t mouse)
{
  view_set_focus(mouse.x, mouse.y, mouse.buttons.right);
}

void onMouseMove(mouse_t mouse)
{
  const float sensitivity = 1.5f;
  static float oldx = 0, oldy = 0;
  if(mouse.buttons.left)
  {
    float x[3] = {1.0f, 0, 0};
    float y[3] = {0, 1.0f, 0};
    float ay = - sensitivity*(mouse.x - oldx)/(float)view_width();
    float ax = sensitivity*(mouse.y - oldy)/(float)view_height();
    view_active_cam_rotate(x, ax);
    view_active_cam_rotate(y, ay);
  }
  oldx = mouse.x;
  oldy = mouse.y;
}

void onClose()
{
  rt.quit = 1;
}


int init(const char *filename, int argc, char *argv[])
{
  // open input file
  FILE *f = fopen(filename, "rb");
  if(!f)
  {
    fprintf(stderr, "[main] init: can't open %s for reading!\n", filename);
    return 1;
  }

  rt.num_threads = sysconf(_SC_NPROCESSORS_ONLN);

  // store cmd line for modules
  rt.argc = argc;
  rt.argv = argv;

  // parse cmd line
  rt.anim_frame = 1;
  rt.batch_frames = 1;      // skip sync point (useful for offline rendering)
  strncpy(rt.output_filename, "render", sizeof(rt.output_filename));
  
  for(int i=0;i<argc;i++)
  {
    if((strcmp(argv[i], "-t") == 0) && (argc > i+1 && argv[i+1][0] != '-')) rt.num_threads = atol(argv[++i]);
    else if((strcmp(argv[i], "-x") == 0)) { if(argc > i+1 && argv[i+1][0] != '-') strncpy(rt.output_filename, argv[++i], sizeof(rt.output_filename)); }
    else if((strcmp(argv[i], "--frame") == 0) && (argc > i+1)) rt.anim_frame = atol(argv[++i]);
    else if((strcmp(argv[i], "--batch") == 0) && (argc > i+1)) rt.batch_frames = atol(argv[++i]);
  }
  fprintf(stdout, "[main] running on %d threads\n", rt.num_threads);

  strncpy(rt.searchpath, filename, 256);
  char *c = rt.searchpath + strlen(rt.searchpath);
  for(;*c != '/' && c != rt.searchpath;c--);
  *c = '\0';
  strncpy(rt.basename, filename, 256);
  c = rt.basename + strlen(rt.basename);
  for(;*c != '.' && c != rt.basename;c--);
  if(c != rt.basename) *c = '\0';


  char datafile[1024];
  // actually this is not what i mean. i mean where the binary lies.
  snprintf(datafile, sizeof(datafile), "data/ergb2spec.coeff");//, rt.searchpath);
  rt.rgb2spec = rgb2spec_init(datafile);
  if(!rt.rgb2spec)
  {
    fprintf(stderr, "[main] could not load `%s', expect trouble!\n", datafile);
    return 1;
  }

  // initialise cameras and resolution, also init rt.display
  rt.view = view_init();
  rt.render = render_init();
  rt.lights = lights_init();

  if(!rt.display)
  {
    fprintf(stderr, "[main] could not open display!\n");
    return 1;
  }

  rt.shader = shader_init(f);

  rt.shader_index = NULL;
  rt.sampler = NULL;
  rt.accel = NULL;

  rt.points = points_init(rt.num_threads, rt.anim_frame);

  rt.prims = (prims_t *)malloc(sizeof(prims_t));
  prims_init(rt.prims);

  rt.threads = threads_init();
  threads_tls_init(rt.threads); // trigger thread local storage initialisation
  rt.sampler = sampler_init();

  rt.play = 0;
  rt.quit = 0;
  rt.screenshot = 0;

  // load scene
  if(common_load_scene(f))
  {
    fprintf(stderr, "[main] could not load nra2 file!\n");
    fclose(f);
    return 1;
  }
  fclose(f);

  // depends on geo, too
  rt.pointsampler = pointsampler_init(rt.anim_frame);

  // build rt accel
  double start, end;
  rt.accel = accel_init(rt.prims);
  printf("[main] creating accel..\n");
  start = common_time_wallclock();
  accel_build(rt.accel, rt.basename);
  end = common_time_wallclock();
  printf("[main] construction took %.3f seconds\n", end - start);

  // scene epsilon (float sucks)
  rt.epsilon = 3e-3f;

  // register gui listeners
  display_register_callbacks(rt.display, &onKeyDown,
      0, &onKeyUp, &onMouseButtonDown, 0, &onMouseMove, 0, 0);

  return 0;
}

void cleanup()
{
  if(rt.display)      display_close(rt.display);
  if(rt.shader)       shader_cleanup(rt.shader);
  if(rt.accel)        accel_cleanup(rt.accel);
  if(rt.shader_index) free(rt.shader_index);
  if(rt.points)       points_cleanup(rt.points);
  if(rt.view)         view_cleanup(rt.view);
  if(rt.render)       render_cleanup(rt.render);
  if(rt.threads)      threads_cleanup(rt.threads);
  if(rt.sampler)      sampler_cleanup(rt.sampler);
  if(rt.pointsampler) pointsampler_cleanup(rt.pointsampler);
  if(rt.prims)        prims_cleanup(rt.prims);
  if(rt.rgb2spec)     rgb2spec_cleanup(rt.rgb2spec);
}

void main_screenshot()
{
  rt.screenshot = 0;

  display_print(rt.display, 0, 0, "saved render");
  printf("[main] saving framebuffers to %s%s\n", rt.basename, rt.output_filename);
  view_write_images(rt.output_filename);
}

void run()
{
  // for timings:
  common_setaffinity(0);

  rt.frames = 0;

  display_print(rt.display, 0, 0, " Tv [12] Av [34] f [56] iso [78] extra controls [tab]");
  const double start_time = common_time_wallclock();
  double start = common_time_wallclock(), end;

  while(!rt.quit)
  {
    view_render();

    display_pump_events(rt.display);
    rt.frames += rt.batch_frames;
    end = common_time_wallclock();
    printf("  %.3f s/frame, %lu spp      \r", end - start, view_overlays());
    start = end;
    fflush(stdout);

    if(rt.screenshot) main_screenshot();
    view_active_cam_move();
  }
  const double end_time = common_time_wallclock();
  fprintf(stdout, "[main] rendered %lu frames in an average of %.3f s/frame\n", rt.frames, (end_time-start_time)/rt.frames);
}

int main(int argc, char* arg[])
{
  char *rev = VERSION;
  printf("corona-13 rev %s\n", rev);
  if(argc == 1)
  {
    fprintf(stderr, "usage: %s <scene.nra2>\n"
                    "  [-c camera[%%04d].cam]                   load camera file, %%04d enables camera motion blur\n"
                    "  [-t threads]\n"
                    "  [-x <postfix>]                          overwrite <basename><postfix>.* screenshot\n"
                    "  [--frame <x>]                           use different random seed for different frame\n"
                    "  [--batch <x>]                           sync framebuffer after x progressions only (improve scaling to many threads)\n"
                    "  [--affinity <file>]                     read specific affinity file (default `affinity')\n", arg[0]);
    view_print_usage();
    display_print_usage();
    exit(1);
  }
  if(init(arg[1], argc-2, arg+2)) exit(2);

  run();
  cleanup();
  exit(0);
}
