/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DISPLAY_GLFW_H
#define DISPLAY_GLFW_H

#include <stdio.h>

// we define that to get the correct mapping in main.c. GLFW passes on key codes, not syms,
// which effectively always assumes an acsii keyboard layout
#define GUI_QWERTZ

// get event declarations
#include "display_common.h"

static inline void display_print_info(FILE *fd)
{
  fprintf(fd, "display  : gl/glfw\n");
}

struct GLFWwindow;
typedef struct display_t
{
	int isShuttingDown;
  int width;
  int height;
  char msg[256];
  int msg_len;
  int msg_x, msg_y;

  unsigned int tex_frame_buffer, program_draw_texture, vao_empty;
  struct GLFWwindow *window;

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

display_t *display_open(const char title[], int width, int height);
int display_update(display_t *d, float* pixels);
void display_close(display_t *d);
void display_pump_events(display_t *d);
void display_print(display_t *d, const int px, const int py, const char* msg, ...);
static inline void display_print_usage() {}
static inline int display_control_add(display_t *d, const char *name, float *storage, float min, float max, float step, int logscale) { return 0; }

#endif
