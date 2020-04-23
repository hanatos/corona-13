#pragma once
#include "rgb2spec.h"
#include "matrix3.h"
#include "colourspaces.h"
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>

#define FRAMEBUFFER_MAGIC 1936686951lu

typedef enum framebuffer_flags_t
{
  FB_XYZ=0,
}
framebuffer_flags_t;

typedef struct framebuffer_header_t
{
  // 8-wide (float buffer after it will be avx aligned)
  uint64_t magic;          // magic number to identify file
  uint64_t width, height;  // dimensions of image
  uint16_t channels;       // floats per pixel
  uint16_t flags;          // identify type of data
  float gain;              // scale factor for e.g. mlt
}
framebuffer_header_t;

typedef struct framebuffer_t
{ // struct to hold runtime pointers, does not end up on disk
  framebuffer_header_t *header;
  float *fb;
  int retain;
  char filename[1024];
}
framebuffer_t;

// map a framebuffer read only
static inline int fb_map(
    framebuffer_t *fb,
    const char *filename)
{
  memset(fb, 0, sizeof(*fb));
  strncpy(fb->filename, filename, sizeof(fb->filename));
  int fd = open(fb->filename, O_RDONLY);
  if(fd == -1)
  {
    snprintf(fb->filename, sizeof(fb->filename), "%s/%s", rt.searchpath, filename);
    fd = open(fb->filename, O_RDONLY);
    if(fd == -1) return 1;
  }

  const uint64_t data_size = lseek(fd, 0, SEEK_END);

  if(data_size < sizeof(framebuffer_header_t)) goto fail;

  fb->header = mmap(0, data_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
  fb->fb = (float *)((uint8_t *)fb->header + sizeof(framebuffer_header_t));

  if(fb->header->magic != FRAMEBUFFER_MAGIC) goto fail;
  if(fb->header->width*fb->header->height*fb->header->channels*sizeof(float) + sizeof(framebuffer_header_t) != data_size) goto fail;
  close(fd);

  fb->retain = 1; // by default, don't delete if we only mapped it, not created it
  return 0;
fail:
  if(fd != -1) close(fd);
  if(fb->header) munmap(fb->header, data_size);
  return 2;
}

// initialise a new framebuffer with file backing
// note that this one will unlink the file once done with it by default.
// change the behaviour by setting fb->retain = 1.
static inline int fb_init(
    framebuffer_t *fb,
    const uint64_t w,
    const uint64_t h,
    const int channels,
    const char *filename)
{
  memset(fb, 0, sizeof(*fb));
  uint64_t size = w * h * channels * sizeof(float);
  int file = open(filename, O_CREAT|O_TRUNC|O_RDWR, 0644);
  ssize_t written = 0;
  while(written < sizeof(framebuffer_header_t))
  {
    ssize_t w = write(file, &fb->header, sizeof(framebuffer_header_t)-written);
    if(w == -1)
    {
      fprintf(stderr, "[framebuffer] ERROR: failed to open frame buffer!\n");
      close(file);
      return 1;
    }
    written += w;
  }
  // seek a bit more and also include header so we can directly map it later
  lseek(file, sizeof(framebuffer_header_t) + size - 1, SEEK_SET);
  written = write(file, &channels, 1);
  if(written != 1) fprintf(stderr, "[framebuffer] ERROR: failed to open frame buffer!\n");
  lseek(file, 0, SEEK_SET);
  fb->header = mmap(0, sizeof(framebuffer_header_t) + size, PROT_READ|PROT_WRITE, MAP_SHARED, file, 0);
  fb->fb = (float *)((uint8_t *)fb->header + sizeof(framebuffer_header_t));
  strncpy(fb->filename, filename, sizeof(fb->filename));
  fb->header->magic = FRAMEBUFFER_MAGIC;
  fb->header->width  = w;
  fb->header->height = h;
  fb->header->channels = channels;
  fb->header->gain = 1.0f;
  close(file);
  return 0;
}

static inline void fb_cleanup(framebuffer_t *fb)
{
  munmap(fb->header, sizeof(framebuffer_header_t) + fb->header->channels * fb->header->width * fb->header->height * sizeof(float));
  if(!fb->retain)
    unlink(fb->filename);
}

static inline void fb_clear(framebuffer_t *fb)
{
  memset(fb->fb, 0, sizeof(float) * fb->header->width*fb->header->height*fb->header->channels);
  fb->header->gain = 1.0f;
}

// make a copy of the file backing with the current state
static inline int fb_save_copy(
    framebuffer_t *fb,
    const char *filename)
{
  size_t data_size = fb->header->width*fb->header->height*fb->header->channels*sizeof(float) + sizeof(framebuffer_header_t);
  int fd = open(filename, O_WRONLY|O_CREAT|O_TRUNC, 0644);
  if(fd < 0) return 1;
  ssize_t w = write(fd, fb->header, data_size);
  close(fd);
  return w != data_size;
}

// will write out a pfm file, accounting for gain
static inline int fb_export(
    framebuffer_t *fb,
    const char *filename,
    const int cbeg,    // first channel to write per pixel
    const int ccnt)    // number of channels to write per pixel (clamped to 3)
{
  FILE* f = fopen(filename, "wb");
  if(!f) return 1;
  // align pfm header to sse, assuming the file will
  // be mmapped to page boundaries.
  char header[1024];
  snprintf(header, 1024, "PF\n%lu %lu\n-1.0", fb->header->width, fb->header->height);
  size_t len = strlen(header);
  fprintf(f, "PF\n%lu %lu\n-1.0", fb->header->width, fb->header->height);
  ssize_t off = 0;
  while((len + 1 + off) & 0xf) off++;
  while(off-- > 0) fprintf(f, "0");
  fprintf(f, "\n");
  for(uint64_t j=0;j<fb->header->height;j++)
  {
    for(uint64_t i=0;i<fb->header->width;i++)
    {
      uint64_t p = fb->header->channels * (i+fb->header->width*j) + cbeg;
      float val[3] = {
        fb->fb[p + (0%ccnt)] * fb->header->gain,
        fb->fb[p + (1%ccnt)] * fb->header->gain,
        fb->fb[p + (2%ccnt)] * fb->header->gain,
      };
      fwrite(val, sizeof(float), 3, f);
    }
  }
  fclose(f);
  return 0;
}

// create a spectral texture framebuffer object from a simple float buffer:
static inline int fb_tex_from_float(
    const rgb2spec_t *r2s,
    const float *data,
    const uint64_t width,
    const uint64_t height,
    const char *fb_file,
    const int hdr)
{
  framebuffer_t fb;
  int stride = hdr ? 4 : 3;
  if(fb_init(&fb, width, height, stride, fb_file))
    return 1;
  double max_mul = 1.0f;
  // for(uint64_t k=0;k<width*height;k++)
  // {
  //   float mul = MAX(MAX(data[3*k+0], data[3*k+1]), data[3*k+2]);
  //   max_mul = MAX(max_mul, mul);
  // }
  for(uint64_t k=0;k<width*height;k++)
  {
    float col[3];
    float mul = MAX(MAX(data[3*k+0], data[3*k+1]), data[3*k+2]);
    for(int i=0;i<3;i++) col[i] = CLAMP(data[3*k+i] / mul, 0.0f, 1.0f);
    rgb2spec_fetch(r2s, col, fb.fb+stride*k);
    if(hdr) fb.fb[stride*k+stride-1] = mul;
  }
  fb.header->gain = max_mul;
  fb.retain = 1;
  fb_cleanup(&fb); // unmap and save data to disk
  return 0;
}

// fetch in integer coordinates on scale (0..w-1, 0..h-1)
static inline const float *fb_fetchi(const framebuffer_t *fb, int i, int j)
{
  if(i < 0 || i >= fb->header->width || j < 0 || j >= fb->header->height) return fb->fb;
  return fb->fb + fb->header->channels*(fb->header->width*j + i);
}

// fetch in normalised texture coordinates in [0,1]^2 (will be repeated if exceeds)
static inline const float *fb_fetch(const framebuffer_t *fb, float s, float t)
{
  const float fu = (s - floorf(s))*(fb->header->width);
  const float fv = (t - floorf(t))*(fb->header->height);
  const int u = (int)fu;
  const int v = (int)fv;
  return fb_fetchi(fb, u, v);
}

// bilinear lookup on [0,1] normalised coordinates
static inline void fb_fetch_bilin(
    const framebuffer_t *fb,
    float s, float t,
    float *res)
{
  const float fu = (s - floorf(s))*fb->header->width;
  const float fv = (t - floorf(t))*fb->header->height;
  const int u = (int)fu;
  const int v = (int)fv;
  const float wu = fu - u;
  const float wv = fv - v;

  // bilinear filter
  const float *v00 = fb_fetchi(fb, u,   v);
  const float *v01 = fb_fetchi(fb, u+1, v);
  const float *v10 = fb_fetchi(fb, u,   v+1);
  const float *v11 = fb_fetchi(fb, u+1, v+1);
  for(int c=0;c<fb->header->channels;c++)
    res[c] = v00[c]*(1-wu)*(1-wv) + v01[c]*wu*(1-wv) + v10[c]*(1-wu)*wv + v11[c]*wu*wv;
}

#undef FRAMEBUFFER_MAGIC
