/*
   This file is part of corona-13.

   copyright (c) 2015 johannes hanika.

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

// fake function:
void view_clear() { }
#include "display.h"
#include "corona_common.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>

#define INBUF_SIZE 4096

// overwrite listener methods.
// buffer for events
#define eventbuf_size 2000
network_event_t eventbuf[eventbuf_size];
int eventbuf_pos = 0;

void onKeyDown(keycode_t code)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EKeyDown;
  eventbuf[eventbuf_pos].data.code = code;
  eventbuf_pos ++;
}

void onKeyPressed(keycode_t code)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EKeyPressed;
  eventbuf[eventbuf_pos].data.code = code;
  eventbuf_pos ++;
}

void onKeyUp(keycode_t code)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EKeyUp;
  eventbuf[eventbuf_pos].data.code = code;
  eventbuf_pos ++;
}

void onMouseButtonUp(mouse_t mouse)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EMouseButtonUp;
  eventbuf[eventbuf_pos].data.mouse = mouse;
  eventbuf_pos ++;
}

void onMouseButtonDown(mouse_t mouse)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EMouseButtonDown;
  eventbuf[eventbuf_pos].data.mouse = mouse;
  eventbuf_pos ++;
}

void onMouseMove(mouse_t mouse)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EMouseMove;
  eventbuf[eventbuf_pos].data.mouse = mouse;
  eventbuf_pos ++;
}
void onActivate(char a)
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EActivate;
  eventbuf[eventbuf_pos].data.activated = a;
  eventbuf_pos ++;
}

void onClose()
{
  if(eventbuf_pos >= eventbuf_size) return;
  eventbuf[eventbuf_pos].type = EClose;
  eventbuf_pos ++;
}


/////////////////////////////////////////////////////////////////
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

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h> 
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

// source functions
void dt_imageio_jpeg_init_source(j_decompress_ptr cinfo) { }
boolean dt_imageio_jpeg_fill_input_buffer(j_decompress_ptr cinfo)
{
  return 1;
}
void dt_imageio_jpeg_skip_input_data(j_decompress_ptr cinfo, long num_bytes)
{
  ssize_t i = cinfo->src->bytes_in_buffer - num_bytes;
  if(i < 0) i = 0;
  cinfo->src->bytes_in_buffer = i;
  cinfo->src->next_input_byte += num_bytes;
}
void dt_imageio_jpeg_term_source(j_decompress_ptr cinfo) { }


int dt_imageio_jpeg_decompress_header(const void *in, size_t length, dt_imageio_jpeg_t *jpg)
{
  jpeg_create_decompress(&(jpg->dinfo));
  jpg->src.init_source = dt_imageio_jpeg_init_source;
  jpg->src.fill_input_buffer = dt_imageio_jpeg_fill_input_buffer;
  jpg->src.skip_input_data = dt_imageio_jpeg_skip_input_data;
  jpg->src.resync_to_restart = jpeg_resync_to_restart;
  jpg->src.term_source = dt_imageio_jpeg_term_source;
  jpg->src.next_input_byte = (JOCTET *)in;
  jpg->src.bytes_in_buffer = length;

  struct dt_imageio_jpeg_error_mgr jerr;
  jpg->dinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = dt_imageio_jpeg_error_exit;
  if(setjmp(jerr.setjmp_buffer))
  {
    jpeg_destroy_decompress(&(jpg->dinfo));
    return 1;
  }

  jpg->dinfo.src = &(jpg->src);
  jpeg_read_header(&(jpg->dinfo), TRUE);
  jpg->dinfo.out_color_space = JCS_RGB;
  jpg->dinfo.out_color_components = 3;
  jpg->width = jpg->dinfo.image_width;
  jpg->height = jpg->dinfo.image_height;
  return 0;
}

static int decompress_plain(dt_imageio_jpeg_t *jpg, uint8_t *out)
{
  JSAMPROW row_pointer[1];
  row_pointer[0] = (uint8_t *)malloc(jpg->dinfo.output_width * jpg->dinfo.num_components);
  uint8_t *tmp = out;
  while(jpg->dinfo.output_scanline < jpg->dinfo.image_height)
  {
    if(jpeg_read_scanlines(&(jpg->dinfo), row_pointer, 1) != 1)
    {
      free(row_pointer[0]);
      return 1;
    }
    for(unsigned int i = 0; i < jpg->dinfo.image_width; i++)
    {
      for(int k = 0; k < 3; k++) tmp[4 * i + k] = row_pointer[0][3 * i + k];
    }
    tmp += 4 * jpg->width;
  }
  free(row_pointer[0]);
  return 0;
}

int dt_imageio_jpeg_decompress(dt_imageio_jpeg_t *jpg, uint8_t *out)
{
  struct dt_imageio_jpeg_error_mgr jerr;
  jpg->dinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = dt_imageio_jpeg_error_exit;
  if(setjmp(jerr.setjmp_buffer))
  {
    jpeg_destroy_decompress(&(jpg->dinfo));
    return 1;
  }

  (void)jpeg_start_decompress(&(jpg->dinfo));

  if(setjmp(jerr.setjmp_buffer))
  {
    jpeg_destroy_decompress(&(jpg->dinfo));
    return 1;
  }

  if(decompress_plain(jpg, out)) return 1;

  if(setjmp(jerr.setjmp_buffer))
  {
    jpeg_destroy_decompress(&(jpg->dinfo));
    return 1;
  }

  (void)jpeg_finish_decompress(&(jpg->dinfo));
  jpeg_destroy_decompress(&(jpg->dinfo));
  return 0;
}
/////////////////////////////////////////////////////////////////


int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "[netrender] usage:  %s <hostname> [port]\n", arg[0]);
    fprintf(stderr, "[netrender] note: if connection fails, you may want to tunnel through by\n");
    fprintf(stderr, "[netrender]   ssh -f hostname -L 8090:hostname:8090 -N\n");
    fprintf(stderr, "[netrender]   %s localhost\n", arg[0]);
    exit(1);
  }
  char *hostname = arg[1];
  int port = 8090;
  if(argc > 2) port = atol(arg[2]);

  int sockfd, udpfd;
  int width, height;
  display_t *display;
  double start, end;

  // open socket
  struct hostent *host = gethostbyname(hostname);
  if (!host)
  {
    fprintf(stderr, "[netrender] unknown host %s!\n", hostname);
    exit(1);
  }

  struct sockaddr_in addr;
  memset(&addr, 0, sizeof(struct sockaddr_in));
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  memcpy(&addr.sin_addr, host->h_addr_list[0], host->h_length);

  printf("[netrender] connecting to %s:%d ...\n", inet_ntoa(addr.sin_addr), port);

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  udpfd  = socket(AF_INET, SOCK_DGRAM, 0);
  if (sockfd < 0 || udpfd < 0)
  {
    perror("socket");
    exit(1);
  }

#ifdef NET_UDP
  // connect to rt server
  struct sockaddr_in serv_addr;
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(port+1);
  if (bind(udpfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0 ||
      connect(sockfd, (struct sockaddr *) &addr, sizeof(struct sockaddr_in)) < 0)
#else
  if (connect(sockfd, (struct sockaddr *) &addr, sizeof(struct sockaddr_in)) < 0)
#endif
  {
    perror("connect");
    exit(1);
  }

  // read width, height.
  {
    char buf[1024];
    read(sockfd, buf, 11);
    if(strncmp(buf, "boom, baby!", 11) != 0)
    {
      fprintf(stderr, "[netrender] got wrong handshake from remote!\n");
      exit(1);
    }
    read(sockfd, &width, sizeof(unsigned int));
    read(sockfd, &height, sizeof(unsigned int));
  }
  printf("[netrender] connected to a %dx%d render stream!\n", width, height);
  uint8_t *jpgbuf = (uint8_t *)malloc(width*height*sizeof(uint32_t));
  uint8_t *framebuf = (uint8_t *)malloc(width*height*sizeof(uint32_t));
  char title[255];
  sprintf(title, "corona on %s", hostname);
  display = display_open(title, width, height);
  display_register_callbacks(display, onKeyDown, onKeyPressed, onKeyUp, onMouseButtonDown, onMouseButtonUp, onMouseMove, onActivate, onClose);

  start = common_time_wallclock();
  int frame = 0;
  while(1)
  {
    int framesize = 0, left = 0;
    uint8_t *inbuf = jpgbuf;
    while(1)
    { // read whole jpg buffer
      int request = framesize ? left : 4096;
#ifdef NET_UDP
      ssize_t size = recv(udpfd, inbuf, request, 0);
#else
      ssize_t size = read(sockfd, inbuf, request);
#endif
      if(size <= 0) goto disconnected;
      if(!framesize)
      {
        framesize = ((uint32_t *)inbuf)[0];
        left = framesize + sizeof(uint32_t);
      }
      left -= size;
      if(left <= 0) break;
      inbuf += size;
    }
    inbuf = jpgbuf + sizeof(uint32_t); // skip bytesize

    dt_imageio_jpeg_t jpg;
    int res = dt_imageio_jpeg_decompress_header(inbuf, framesize, &jpg);
    res = dt_imageio_jpeg_decompress(&jpg, framebuf);

    display_update_rgba(display, (unsigned int*)framebuf);
    display_pump_events(display);
    // send events
    if(eventbuf_pos)
    {
      write(sockfd, eventbuf, eventbuf_pos * sizeof(network_event_t));
      eventbuf_pos = 0;
    }

    if((frame++ & 7) == 7)
    {
      end = common_time_wallclock();
      printf(" %0.2f s/frame    \r", (end - start)/8.0);
      fflush(stdout);
      start = end;
    }
  }
disconnected:
  close(sockfd);
  close(udpfd);
  free(jpgbuf);
  free(framebuf);
  display_close(display);
  printf("[netrender] disconnected.\n");

  return 0;
}
