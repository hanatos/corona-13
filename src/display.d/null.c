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
#include "display.h"
#include "view.h"

#include <stdarg.h>

void display_print_info(FILE *fd)
{
  fprintf(fd, "display  : null (render to file)\n");
}

typedef struct display_t
{
  uint64_t frames, backup_frames, sequence_frames;
  double timeout, start_time;
  void (*onKeyDown)(keycode_t);
  void (*onKeyPressed)(keycode_t);
  void (*onKeyUp)(keycode_t);
  void (*onMouseButtonUp)(mouse_t);
  void (*onMouseButtonDown)(mouse_t);
  void (*onMouseMove)(mouse_t);
  void (*onActivate)(char);
  void (*onClose)();
}
display_t;

display_t *display_open(const char title[], int width, int height)
{
  display_t *p = (display_t *)malloc(sizeof(display_t));
  p->frames = 10;
  p->timeout = 0.0;
  p->backup_frames = 0;
  p->sequence_frames = 0;
  for(int k=0;k<rt.argc;k++) if(rt.argv[k][0] == '-' && rt.argv[k][1] == 's' && k < rt.argc-1) p->frames = atol(rt.argv[++k]);
                        else if(rt.argv[k][0] == '-' && rt.argv[k][1] == 'b' && k < rt.argc-1) p->backup_frames = atol(rt.argv[++k]);
                        else if(rt.argv[k][0] == '-' && rt.argv[k][1] == 'q') p->sequence_frames = 1;
                        else if(rt.argv[k][0] == '-' && rt.argv[k][1] == 'o' && k < rt.argc-1) p->timeout = atof(rt.argv[++k]);
  printf("[display] simulating %lu samples per pixel, backup image every %lu spp,\n"
         "          timeout after %.2f seconds,\n",
         p->frames, p->backup_frames, p->timeout);
  if(p->sequence_frames)
    printf("          sequence shot every power-of-two progressions.\n");
  p->start_time = common_time_wallclock();
  return p;
}
int display_update(display_t *d, const float *pixels, const float scale)
{
  if(d->backup_frames && (((rt.frames+rt.batch_frames) % d->backup_frames) == 0)) view_write_images("_backup");
  // this checks if it's a power of 2
  const uint64_t n = rt.frames+rt.batch_frames;
  if(d->sequence_frames && !(n == 0 || (n & (n-1))))
  {
    char filename[1024];
    snprintf(filename, sizeof(filename), "%s_seq_%06lu", rt.output_filename, n);
    view_write_images(filename);
  }
  if((d->timeout == 0.0f) && rt.frames + rt.batch_frames >= d->frames) rt.screenshot = rt.quit = 1;
  if((d->timeout > 0.0f) && (common_time_wallclock() - d->start_time > d->timeout)) { d->timeout = 0.0; view_write_images(rt.output_filename); rt.quit = 1; }
  return 0;
}
int display_update_rgba(display_t *d, const unsigned int* rgba) { return 0; }
void display_close(display_t *d) {}
void display_pump_events(display_t *d) {}
void display_print(display_t *d, const int px, const int py, const char* msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf("[display] ");
  vprintf(msg, ap);
  printf("\n");
  va_end(ap);
}
void display_print_usage()
{
  fprintf(stderr, "  [-s samples]                            write image after number of samples and quit\n"
                  "  [-b backup_frames]                      write single backup frame every backup_frames progressions\n"
                  "  [-q ]                                   write and keep sequence of frames every power-of-two progressions\n"
                  "  [-o timeout]                            write image after timeout seconds and quit (overrides -s)\n");
}
int display_control_add(display_t *d, const char *name, float *storage, float min, float max, float step, int logscale, int clear) { return 1; }

void display_register_callbacks(display_t *d, 
  void (*onKeyDown)(keycode_t),
  void (*onKeyPressed)(keycode_t),
  void (*onKeyUp)(keycode_t),
  void (*onMouseButtonDown)(mouse_t),
  void (*onMouseButtonUp)(mouse_t),
  void (*onMouseMove)(mouse_t),
  void (*onActivate)(char),
  void (*onClose)())
{}
int display_requires_update(display_t *d) { return 0; }
