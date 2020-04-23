#include "corona_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>

#include <sys/stat.h>
#include <sys/mman.h>

#define MAX_DBORS 20

typedef struct dbor_buffer_s
{
  int fd;              // file descriptor
  void *data;          // mapped buffer
  size_t data_size;    // size of file
  
  int width;
  int height;
  float* pixel;
}
dbor_buffer_t;


int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "[dbor] usage: input_file [K_min] [K].\n");
    fprintf(stderr, "       K_min: minimum number of samples/pixel, default 0.01\n");
    fprintf(stderr, "       K    : used for reweighting, set to -large to avoid reweighting. default 1\n");
    exit(1);
  }

  float K_min = argc > 2 ? atof(arg[2]) : 0.01;
  float K = argc > 3 ? atof(arg[3]) : 10;
    
  // load dbor cascade
  int num_dbors = 0;
  dbor_buffer_t dbor_cascade[MAX_DBORS] = {{ 0 }};
  while (1)
  {
    if(num_dbors >= MAX_DBORS)
    {
      fprintf(stderr, "[dbor] MAX_DBORS of %d reached.\n", MAX_DBORS);
      break;
    }

    char filename[1024];
    snprintf(filename, sizeof(char)*1024, "%s_dbor%02d.pfm", arg[1], num_dbors);
    if(access(filename, F_OK))
    {
      break;
    }
    else
    {
      dbor_buffer_t* db = &dbor_cascade[num_dbors];

      db->fd = open(filename, O_RDONLY);
      if (db->fd == -1) {
        fprintf(stderr, "[dbor] could not open dbor cascade buffer `%s'\n", filename);
        continue;
      }

      // mmap data
      db->data_size = lseek(db->fd, 0, SEEK_END);
      lseek(db->fd, 0, SEEK_SET);
      db->data = mmap(0, db->data_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, db->fd, 0);
      
      // get pfm header for faster grabbing later on.
      // no error handling is done, no comments supported.
      char *endptr;
      db->width = strtol(db->data+3, &endptr, 10);
      db->height = strtol(endptr, &endptr, 10);
      endptr++; // remove newline
      while(endptr < (char *)db->data + 100 && *endptr != '\n') endptr++;
      db->pixel = (float*)(++endptr); // remove second newline

      // get pfm meta data
      // printf("file size %ld, width %d, height %d\n", db->data_size, db->width, db->height);

      num_dbors++;
    }
  }
  // printf("num dbors %d\n", num_dbors);

  // create file for merged dbor
  dbor_buffer_t dm = { 0 };
  char filename[1024];
  snprintf(filename, sizeof(filename), "%s_dbor.pfm", arg[1]);
  dm.fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0644);
  if (dm.fd == -1) {
    fprintf(stderr, "[dbor] could not open dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  dm.data_size = dbor_cascade[0].data_size;
  if (lseek(dm.fd, dm.data_size-1, SEEK_SET) == -1) {
    fprintf(stderr, "[dbor] could not reserve memory for dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  if (write(dm.fd, "", 1) == -1) {
    fprintf(stderr, "[dbor] could write memory for dbor merged buffer `%s'\n", filename);
    goto fail;
  }
  lseek(dm.fd, 0, SEEK_SET);
  dm.data = mmap(0, dm.data_size, PROT_WRITE | PROT_READ, MAP_SHARED | MAP_NORESERVE, dm.fd, 0);
  memcpy(dm.data, dbor_cascade[1].data, dm.data_size);
  fprintf(stdout, "[dbor] writing `%s'\n", filename);
  
  // get pfm header for checking memcpy
  char *endptr;
  dm.width = strtol(dm.data+3, &endptr, 10);
  dm.height = strtol(endptr, &endptr, 10);
  endptr++; // remove newline
  while(endptr < (char *)dm.data + 100 && *endptr != '\n') endptr++;
  dm.pixel = (float*)(++endptr); // remove second newline
  assert(dm.width == dbor_cascade[0].width);
  assert(dm.height == dbor_cascade[0].height);
  // printf("result: file size %ld, width %d, height %d\n", dm.data_size, dm.width, dm.height);

  for (int h = 0; h < dm.height; ++h)
  for (int w = 0; w < dm.width; ++w)
  {
    const int offset = 3*(h*dm.width+w);
    
    for (int c = 0; c < 3; ++c) 
      dm.pixel[offset+c] = dbor_cascade[0].pixel[offset+c];
#if 0
    for (int i = 1; i < num_dbors; ++i)
        for (int c = 0; c < 3; ++c)
          dm.pixel[offset+c] += dbor_cascade[i].pixel[offset+c];
#else
    for (int i = 1; i < num_dbors; ++i)
    {
      float n = 0, n_avg = 0;
      for (int dh = -1; dh < 2; ++dh)
      for (int dw = -1; dw < 2; ++dw)
      {
        const int offset = 3*(CLAMP(h+dh, 0, dm.height-1) * dm.width + CLAMP(w+dw, 0, dm.width-1));
        float n_tmp = 0;
        if (i == (num_dbors - 1)) 
        {
          float* p0 = dbor_cascade[i-1].pixel;
          float* p1 = dbor_cascade[i].pixel;
          n_tmp = 1.f/3.f*(p0[offset+0]+p0[offset+1]+p0[offset+2]) / (1 << (i-1))
                + 1.f/3.f*(p1[offset+0]+p1[offset+1]+p1[offset+2]) / (1 << i);
        }
        else {
          float* p0 = dbor_cascade[i-1].pixel;
          float* p1 = dbor_cascade[i].pixel;
          float* p2 = dbor_cascade[i+1].pixel;
          n_tmp = 1.f/3.f*(p0[offset+0]+p0[offset+1]+p0[offset+2]) / (1 << (i-1))
                + 1.f/3.f*(p1[offset+0]+p1[offset+1]+p1[offset+2]) / (1 << i)
                + 1.f/3.f*(p2[offset+0]+p2[offset+1]+p2[offset+2]) / (1 << (i+1));
        }
        
        n_avg += n_tmp;
        if (dh == 0 && dw == 0)
          n = n_tmp;
      }
      n_avg /= 9.f;
      
      if (n_avg > K_min && n > K_min)
      {
        const float weight = n < (K + K_min) ? (n - K_min)/K : 1.f;
        for (int c = 0; c < 3; ++c)
          dm.pixel[offset+c] += weight * dbor_cascade[i].pixel[offset+c];
      }
    }
#endif
  }

fail:
  for (int i = 0; i < num_dbors; ++i)
  {
    dbor_buffer_t* db = &dbor_cascade[i];
    munmap(db->data, db->data_size);
    close(db->fd);
  }
  munmap(dm.data, dm.data_size);
  close(dm.fd);
}
