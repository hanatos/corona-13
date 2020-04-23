#pragma once
#include "prims.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

// super stupid sidecar data with an array of floats per vertex.
// could encode anything in it.

// TODO: some more clever knowledge about compressed types/channels per vertex?
typedef struct vdata_t
{
  int num_verts;
  int num_slots;
  int fd;
  float *data;
  size_t data_size;
}
vdata_t;

static inline int vdata_map(vdata_t *d, const char *filename)
{
  memset(d, 0, sizeof(*d));
  d->fd = open(filename, O_RDONLY);
  if(d->fd == -1) return 1;
  d->data_size = lseek(d->fd, 0, SEEK_END);
  lseek(d->fd, 0, SEEK_SET);
  d->data = mmap(0, d->data_size, PROT_READ, MAP_SHARED, d->fd, 0);
  d->num_verts = d->data_size / sizeof(float);
  d->num_slots = 1; // <= obviously a lie
  return 0;
}

static inline void vdata_alloc(vdata_t *d, const int num_verts, const int num_slots)
{
  memset(d, 0, sizeof(*d));
  d->fd = -1;
  d->data_size = num_verts*num_slots*sizeof(float);
  d->num_verts = num_verts;
  d->num_slots = num_slots;
  d->data = aligned_alloc(32, d->data_size);
}

static inline int vdata_write(const vdata_t *d, const char *filename)
{
  assert(d->fd < 0);
  int fd = open(filename, O_WRONLY|O_CREAT|O_TRUNC, 0644);
  if(fd < 0) return 1;
  ssize_t w = write(fd, d->data, d->data_size);
  close(fd);
  return w != d->data_size;
}

static inline void vdata_cleanup(vdata_t *d)
{
  if(d->fd >= 0)
  {
    munmap(d->data, d->data_size);
    close(d->fd);
  }
  else free(d->data);
  memset(d, 0, sizeof(*d));
}

static inline float *vdata_get(vdata_t *d, int v)
{
  assert((d->num_slots * v)/sizeof(float) < d->data_size);
  return d->data + d->num_slots * v;
}
