/*
   This file is part of corona-13.

   copyright (c) 2015-2016 johannes hanika.

   corona-13 is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   corona-13 is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "corona_common.h"
#include "display.h"
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <sys/poll.h>
#include <strings.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <assert.h>

void display_print_info(FILE *fd)
{
  fprintf(fd, "display  : network (mjpeg stream)\n");
}

typedef struct display_t
{
  int sock;
  int acceptsock;
  int udpsock;
  unsigned int *buffer, *sendbuf;
  int width;
  int height;

  char msg[256];
  int msg_len;
  int msg_x, msg_y;

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

void display_print_usage()
{
  fprintf(stderr, "  [-p port]                               bind to different port number (default 8090)\n");
}

int display_control_add(display_t *d, const char *name, float *storage, float min, float max, float step, int logscale, int clear) { return 1; }


network_event_t event;
size_t remaining;

double resolution;

int keyMapsInitialized = 0;
const unsigned char font9x16[] = {
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,30,0,63,176,63,176,30,0,0,0,0,0,0,0,0,0,112,0,120,0,0,0,0,0,120,0,112,0,0,0,0,0,4,64,31,240,31,240,4,64,
4,64,31,240,31,240,4,64,0,0,28,96,62,48,34,16,226,28,226,28,51,240,25,224,0,0,0,0,24,48,24,96,0,192,1,128,3,0,6,0,12,48,24,48,0,0,1,224,27,240,62,16,39,16,61,224,27,240,2,16,0,0,0,0,
0,0,8,0,120,0,112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,192,31,224,48,48,32,16,0,0,0,0,0,0,0,0,0,0,32,16,48,48,31,224,15,192,0,0,0,0,0,0,1,0,5,64,7,192,3,128,3,128,
7,192,5,64,1,0,0,0,0,0,1,0,1,0,7,192,7,192,1,0,1,0,0,0,0,0,0,0,0,0,0,8,0,120,0,112,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,48,0,48,0,0,0,0,0,0,0,0,0,48,0,96,0,192,1,128,3,0,6,0,12,0,0,0,0,0,15,192,31,224,48,48,35,16,48,48,31,224,15,192,0,0,0,0,0,0,8,16,24,16,63,240,63,240,0,16,
0,16,0,0,0,0,16,112,48,240,33,144,35,16,38,16,60,48,24,48,0,0,0,0,16,32,48,48,32,16,34,16,34,16,63,240,29,224,0,0,0,0,3,0,7,0,13,0,25,16,63,240,63,240,1,16,0,0,0,0,62,32,62,48,
34,16,34,16,34,16,35,240,33,224,0,0,0,0,15,224,31,240,50,16,34,16,34,16,35,240,1,224,0,0,0,0,48,0,48,0,33,240,35,240,38,0,60,0,56,0,0,0,0,0,29,224,63,240,34,16,34,16,34,16,63,240,29,224,
0,0,0,0,28,0,62,16,34,16,34,16,34,48,63,224,31,192,0,0,0,0,0,0,0,0,0,0,12,96,12,96,0,0,0,0,0,0,0,0,0,0,0,0,0,16,12,112,12,96,0,0,0,0,0,0,0,0,0,0,1,0,3,128,
6,192,12,96,24,48,16,16,0,0,0,0,0,0,4,128,4,128,4,128,4,128,4,128,4,128,0,0,0,0,0,0,16,16,24,48,12,96,6,192,3,128,1,0,0,0,0,0,24,0,56,0,32,0,33,176,35,176,62,0,28,0,0,0,
0,0,15,224,31,240,16,16,19,208,19,208,31,208,15,128,0,0,0,0,7,240,15,240,25,0,49,0,25,0,15,240,7,240,0,0,0,0,32,16,63,240,63,240,34,16,34,16,63,240,29,224,0,0,0,0,15,192,31,224,48,48,32,16,
32,16,48,48,24,96,0,0,0,0,32,16,63,240,63,240,32,16,48,48,31,224,15,192,0,0,0,0,32,16,63,240,63,240,34,16,39,16,48,48,56,112,0,0,0,0,32,16,63,240,63,240,34,16,39,0,48,0,56,0,0,0,0,0,
15,192,31,224,48,48,33,16,33,16,49,224,25,240,0,0,0,0,63,240,63,240,2,0,2,0,2,0,63,240,63,240,0,0,0,0,0,0,0,0,32,16,63,240,63,240,32,16,0,0,0,0,0,0,0,224,0,240,0,16,32,16,63,240,
63,224,32,0,0,0,0,0,32,16,63,240,63,240,3,0,7,128,60,240,56,112,0,0,0,0,32,16,63,240,63,240,32,16,0,16,0,48,0,112,0,0,0,0,63,240,63,240,28,0,14,0,28,0,63,240,63,240,0,0,0,0,63,240,
63,240,28,0,14,0,7,0,63,240,63,240,0,0,0,0,31,224,63,240,32,16,32,16,32,16,63,240,31,224,0,0,0,0,32,16,63,240,63,240,34,16,34,0,62,0,28,0,0,0,0,0,31,224,63,240,32,16,32,48,32,28,63,252,
31,228,0,0,0,0,32,16,63,240,63,240,34,0,35,0,63,240,28,240,0,0,0,0,24,96,60,112,38,16,34,16,35,16,57,240,24,224,0,0,0,0,0,0,56,0,48,16,63,240,63,240,48,16,56,0,0,0,0,0,63,224,63,240,
0,16,0,16,0,16,63,240,63,224,0,0,0,0,63,128,63,192,0,96,0,48,0,96,63,192,63,128,0,0,0,0,63,224,63,240,0,112,3,192,0,112,63,240,63,224,0,0,0,0,48,48,60,240,15,192,7,128,15,192,60,240,48,48,
0,0,0,0,0,0,60,0,62,16,3,240,3,240,62,16,60,0,0,0,0,0,56,112,48,240,33,144,35,16,38,16,60,48,56,112,0,0,0,0,0,0,0,0,63,240,63,240,32,16,32,16,0,0,0,0,0,0,28,0,14,0,7,0,
3,128,1,192,0,224,0,112,0,0,0,0,0,0,0,0,32,16,32,16,63,240,63,240,0,0,0,0,0,0,16,0,48,0,96,0,192,0,96,0,48,0,16,0,0,0,0,0,0,4,0,4,0,4,0,4,0,4,0,4,0,4,0,4,
0,0,0,0,0,0,96,0,112,0,16,0,0,0,0,0,0,0,0,0,0,224,5,240,5,16,5,16,7,224,3,240,0,16,0,0,0,0,32,16,63,240,63,224,4,16,6,16,3,240,1,224,0,0,0,0,3,224,7,240,4,16,4,16,
4,16,6,48,2,32,0,0,0,0,1,224,3,240,6,16,36,16,63,224,63,240,0,16,0,0,0,0,3,224,7,240,5,16,5,16,5,16,7,48,3,32,0,0,0,0,0,0,2,16,31,240,63,240,34,16,48,0,24,0,0,0,0,0,
3,228,7,246,4,18,4,18,3,254,7,252,4,0,0,0,0,0,32,16,63,240,63,240,2,0,4,0,7,240,3,240,0,0,0,0,0,0,0,0,4,16,55,240,55,240,0,16,0,0,0,0,0,0,0,0,0,4,0,6,0,2,4,2,
55,254,55,252,0,0,0,0,32,16,63,240,63,240,1,128,3,192,6,112,4,48,0,0,0,0,0,0,0,0,32,16,63,240,63,240,0,16,0,0,0,0,0,0,7,240,7,240,6,0,3,240,3,240,6,0,7,240,3,240,0,0,4,0,
7,240,3,240,4,0,4,0,7,240,3,240,0,0,0,0,3,224,7,240,4,16,4,16,4,16,7,240,3,224,0,0,0,0,4,2,7,254,3,254,4,18,4,16,7,240,3,224,0,0,0,0,3,224,7,240,4,16,4,18,3,254,7,254,
4,2,0,0,0,0,4,16,7,240,3,240,6,16,4,0,6,0,2,0,0,0,0,0,3,32,7,176,4,144,4,144,4,144,6,240,2,96,0,0,0,0,4,0,4,0,31,224,63,240,4,16,4,48,0,32,0,0,0,0,7,224,7,240,
0,16,0,16,7,224,7,240,0,16,0,0,0,0,7,192,7,224,0,48,0,16,0,48,7,224,7,192,0,0,0,0,7,224,7,240,0,48,0,224,0,224,0,48,7,240,7,224,0,0,4,16,6,48,3,96,1,192,1,192,3,96,6,48,
4,16,0,0,7,226,7,242,0,18,0,18,0,22,7,252,7,248,0,0,0,0,6,48,6,112,4,208,5,144,7,16,6,48,4,48,0,0,0,0,0,0,2,0,2,0,31,224,61,240,32,16,32,16,0,0,0,0,0,0,0,0,0,0,
62,248,62,248,0,0,0,0,0,0,0,0,0,0,32,16,32,16,61,240,31,224,2,0,2,0,0,0,0,0,32,0,96,0,64,0,96,0,32,0,96,0,64,0,0,0,0,0,1,224,3,224,6,32,12,32,6,32,3,224,1,224,0,0};

display_t *display_open(const char title[], int width, int height)
{
  // no x display, just avcodec output

  int port = 8090;
  for(int k=0;k<rt.argc;k++) if(rt.argv[k][0] == '-' && rt.argv[k][1] == 'p' && k < rt.argc-1) port = atol(rt.argv[++k]);

  display_t *d = (display_t *)malloc(sizeof(display_t));
  d->onKeyUp = NULL;
  d->onKeyDown = NULL;
  d->onKeyPressed = NULL;
  d->onMouseButtonDown = NULL;
  d->onMouseButtonUp = NULL;
  d->onMouseMove = NULL;
  d->onActivate = NULL;
  d->onClose = NULL;
  d->width = width;
  d->height = height;
  d->udpsock = -1;

  socklen_t clilen;
  struct sockaddr_in serv_addr, cli_addr;
  d->acceptsock = socket(AF_INET, SOCK_STREAM, 0);
  //d->acceptsock = socket(AF_INET, SOCK_DGRAM, 0);
  if (d->acceptsock < 0)
  {
    perror("opening socket");
    return NULL;
  }
  int opt = 1;
  setsockopt(d->acceptsock, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
#if 1
  memset((char *) &serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = htons(INADDR_ANY);
  serv_addr.sin_port = htons(port);
  if (bind(d->acceptsock, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
  {
    perror("binding");
    return NULL;
  }
  listen(d->acceptsock, 5);
  clilen = sizeof(cli_addr);
  printf("[netrender] waiting for corona-netrender to connect on port %d...\n", port);
  d->sock = accept(d->acceptsock, (struct sockaddr *) &cli_addr, &clilen);
  if (d->sock < 0)
  {
    perror("accept");
    return NULL;
  }
#else
  d->sock = d->acceptsock;
#endif

  char buf[1024];
  sprintf(buf, "boom, baby!");
  int ret = 0;
  ret += write(d->sock, buf, 11);
  ret += write(d->sock, &width, sizeof(unsigned int));
  ret += write(d->sock, &height, sizeof(unsigned int));
  assert(ret == ret);

#ifdef NET_UDP
  // start udp stream
  d->udpsock = socket(AF_INET, SOCK_DGRAM, 0);
  if (d->udpsock < 0)
  {
    perror("opening socket");
    return NULL;
  }
  cli_addr.sin_port = htons(port+1);
  if (connect(d->udpsock, (struct sockaddr *) &cli_addr, clilen) < 0)
  {
    perror("connect");
    return NULL;
  }
#endif

  d->buffer = (unsigned int*) common_alloc(128, sizeof(unsigned int) * width * height);
  d->sendbuf = (unsigned int*) common_alloc(128, sizeof(unsigned int) * width * height);

  return d;
}

void display_close(display_t *d)
{	
  close(d->sock);
  close(d->acceptsock);
  if(d->udpsock >= 0) close(d->udpsock);

  free(d->buffer);
  free(d->sendbuf);
  free(d);
}

void display_pump_events(display_t *d)
{
  // parse event from event sock, pass to listener methods
  struct pollfd fds;
  fds.fd = d->sock;
  fds.events = POLLIN;
  fds.revents = 0;

  while(poll(&fds, 1, 4) > 0 && (fds.revents & POLLIN))
  {
    int size = read(d->sock, &event + sizeof(network_event_t) - remaining, remaining);
    if(size < 0) return;
    remaining -= size;
    if(remaining > 0) return;
    remaining = sizeof(network_event_t);
    switch(event.type)
    {
      case EKeyDown:
        if(d->onKeyDown) d->onKeyDown(event.data.code);
        break;
      case EKeyPressed:
        if(d->onKeyPressed) d->onKeyPressed(event.data.code);
        break;
      case EKeyUp:
        if(d->onKeyUp) d->onKeyUp(event.data.code);
        break;
      case EMouseButtonUp:
        if(d->onMouseButtonUp) d->onMouseButtonUp(event.data.mouse);
        break;
      case EMouseButtonDown:
        if(d->onMouseButtonDown) d->onMouseButtonDown(event.data.mouse);
        break;
      case EMouseMove:
        if(d->onMouseMove) d->onMouseMove(event.data.mouse);
        break;
      case EActivate:
        if(d->onActivate) d->onActivate(event.data.activated);
        break;
      case EClose:
        if(d->onClose) d->onClose();
        break;
      default:
        fprintf(stderr, "[netrender] unknown event %d!\n", event.type);
        break;
    }
  }
}




/////////////////////////////////////////////////////////////////////////////////////

/*
    This file is part of darktable,
    copyright (c) 2009--2010 johannes hanika.
    copyright (c) 2015 LebedevRI.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
    (stolen from darktable and heavily cannibalised for corona)
*/

// this fixes a rather annoying, long time bug in libjpeg :(
#undef HAVE_STDLIB_H
#undef HAVE_STDDEF_H
#include <jpeglib.h>
#undef HAVE_STDLIB_H
#undef HAVE_STDDEF_H
#include <setjmp.h>
  
typedef struct dt_imageio_jpeg_t
{ 
  int width, height;
  struct jpeg_source_mgr src;
  struct jpeg_destination_mgr dest;
  struct jpeg_decompress_struct dinfo;
  struct jpeg_compress_struct cinfo;
}
dt_imageio_jpeg_t;


// error functions

struct dt_imageio_jpeg_error_mgr
{
  struct jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
}
dt_imageio_jpeg_error_mgr;

typedef struct dt_imageio_jpeg_error_mgr *dt_imageio_jpeg_error_ptr;

void dt_imageio_jpeg_error_exit(j_common_ptr cinfo)
{
  dt_imageio_jpeg_error_ptr myerr = (dt_imageio_jpeg_error_ptr)cinfo->err;
  (*cinfo->err->output_message)(cinfo);
  longjmp(myerr->setjmp_buffer, 1);
}

// destination functions
void dt_imageio_jpeg_init_destination(j_compress_ptr cinfo) { }
boolean dt_imageio_jpeg_empty_output_buffer(j_compress_ptr cinfo)
{
  fprintf(stderr, "[imageio_jpeg] output buffer full!\n");
  return FALSE;
}
void dt_imageio_jpeg_term_destination(j_compress_ptr cinfo) { }

int dt_imageio_jpeg_compress(const uint8_t *in, uint8_t *out, const int width, const int height,
                             const int quality)
{
  struct dt_imageio_jpeg_error_mgr jerr;
  dt_imageio_jpeg_t jpg;
  jpg.dest.init_destination = dt_imageio_jpeg_init_destination;
  jpg.dest.empty_output_buffer = dt_imageio_jpeg_empty_output_buffer;
  jpg.dest.term_destination = dt_imageio_jpeg_term_destination;
  jpg.dest.next_output_byte = (JOCTET *)out;
  jpg.dest.free_in_buffer = 4 * width * height * sizeof(uint8_t);

  jpg.cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = dt_imageio_jpeg_error_exit;
  if(setjmp(jerr.setjmp_buffer))
  {
    jpeg_destroy_compress(&(jpg.cinfo));
    return 1;
  }
  jpeg_create_compress(&(jpg.cinfo));
  jpg.cinfo.dest = &(jpg.dest);

  jpg.cinfo.image_width = width;
  jpg.cinfo.image_height = height;
  jpg.cinfo.input_components = 3;
  jpg.cinfo.in_color_space = JCS_RGB;
  jpeg_set_defaults(&(jpg.cinfo));
  jpeg_set_quality(&(jpg.cinfo), quality, TRUE);
  if(quality > 90) jpg.cinfo.comp_info[0].v_samp_factor = 1;
  if(quality > 92) jpg.cinfo.comp_info[0].h_samp_factor = 1;
  jpeg_start_compress(&(jpg.cinfo), TRUE);
  uint8_t row[3 * width];
  const uint8_t *buf;
  while(jpg.cinfo.next_scanline < jpg.cinfo.image_height)
  {
    JSAMPROW tmp[1];
    buf = in + jpg.cinfo.next_scanline * jpg.cinfo.image_width * 4;
    for(int i = 0; i < width; i++)
      for(int k = 0; k < 3; k++) row[3 * i + k] = buf[4 * i + k];
    tmp[0] = row;
    jpeg_write_scanlines(&(jpg.cinfo), tmp, 1);
  }
  jpeg_finish_compress(&(jpg.cinfo));
  jpeg_destroy_compress(&(jpg.cinfo));
  return 4 * width * height * sizeof(uint8_t) - jpg.dest.free_in_buffer;
}

/////////////////////////////////////////////////////////////////////////////////////


int display_update(display_t *d, const float *pixels)
{
  // encode frame.
  if (rt.quit)
  {
    display_close(d);
    return 0;
  }

  const int w = d->width;
  const int h = d->height;
  for ( unsigned int i = 0; i < w*h; ++i)
  {
    const uint32_t r = fminf(255.999f, fmaxf(0.0f, 256*pixels[4*i+1]));
    const uint32_t g = fminf(255.999f, fmaxf(0.0f, 256*pixels[4*i+2]));
    const uint32_t b = fminf(255.999f, fmaxf(0.0f, 256*pixels[4*i+3]));
    // argb
    d->buffer[i] = (r<<16)| (g<<8) | b;
  }
  // render message:
  int px = d->msg_x;
  for (int pos = 0; pos < d->msg_len; pos++)
  {
    int charPos = (d->msg[pos] - 32)*9*2;
    for (int x = 0; x < 9*2; x++)
    {
      if (px >= d->width) goto render_message_out;
      unsigned char cLine = font9x16[charPos+x];			
      for (int i = 0; i < 8; i++)
      {
        int y = d->msg_y - (15 - (i + (x&1)*8));
        if ((y >= d->height) || (y < 0)) goto render_message_out;
        if (cLine & (1<<(7-i)))
          for(int k=1;k<4;k++) ((unsigned char *)d->buffer)[(px + y*d->width)*4+k] |= 0x80;
        else
          for(int k=1;k<4;k++) ((unsigned char *)d->buffer)[(px + y*d->width)*4+k] >>= 1;
      }
      if (x&1) px++;
    }
  }
render_message_out:;

  int size = dt_imageio_jpeg_compress((uint8_t*)d->buffer, ((uint8_t*)d->sendbuf) + 4, w, h, 90);
  if(size <= 0) return 1;
  // fprintf(stderr, "compressed frame size %d\n", size);
  d->sendbuf[0] = size;
  size += 4; // bufsize

#ifdef NET_UDP
  const int packet_size = 4096;
  for(int i=0;i<size;i+=packet_size)
  {
    int rest = packet_size;
    if(size - i < packet_size) rest = size - i;
    send(d->udpsock, d->sendbuf + i, rest, 0);
  }
  return 1;
#else
  int ret = write(d->sock, d->sendbuf, size);
  return ret == -1;
#endif
}

void display_print(display_t *d, const int px, const int py, const char *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  vsnprintf(d->msg, 255, msg, ap);
  // vprintf(msg, ap);
  // printf("\n");
  va_end(ap);
  d->msg_len = strlen(d->msg);
  d->msg_x = px;
  d->msg_y = py + 15;
}

void display_register_callbacks(display_t *d, 
  void (*onKeyDown)(keycode_t),
  void (*onKeyPressed)(keycode_t),
  void (*onKeyUp)(keycode_t),
  void (*onMouseButtonDown)(mouse_t),
  void (*onMouseButtonUp)(mouse_t),
  void (*onMouseMove)(mouse_t),
  void (*onActivate)(char),
  void (*onClose)())
{
  d->onKeyDown = onKeyDown;
  d->onKeyPressed = onKeyPressed;
  d->onKeyUp = onKeyUp;
  d->onMouseButtonDown = onMouseButtonDown;
  d->onMouseButtonUp = onMouseButtonUp;
  d->onMouseMove = onMouseMove;
  d->onActivate = onActivate;
  d->onClose = onClose;
}

int display_requires_update(display_t *d) { return 1; }
