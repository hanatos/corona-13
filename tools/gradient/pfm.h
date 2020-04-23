#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct pfm_t
{
  float *pixel;
  int64_t wd, ht;
}
pfm_t;

static inline int pfm_load(pfm_t *img, const char *filename)
{ 
  FILE *f = fopen(filename, "rb");
  if(!f) return 1;
  if(fscanf(f, "PF\n%ld %ld\n%*[^\n]", &img->wd, &img->ht) != 2) return 2;
  (void)fgetc(f);
  img->pixel = (float *)malloc(sizeof(float)*3*img->wd*img->ht);
  fread(img->pixel, sizeof(float)*3, img->wd*img->ht, f);
  for(int k=0;k<img->wd*img->ht*3;k++)
    // if(!(img->pixel[k]==img->pixel[k])) img->pixel[k] = 0;
    if(!(img->pixel[k]< 10e6f)) img->pixel[k] = 0;
  fclose(f);
  return 0;
} 

static inline int pfm_write(const pfm_t *img, const char *filename)
{ 
  FILE *f = fopen(filename, "wb");
  if(!f) return 1;
  char header[1024];
  snprintf(header, 1024, "PF\n%lu %lu\n-1.0", img->wd, img->ht);
  size_t len = strlen(header);
  fprintf(f, "PF\n%lu %lu\n-1.0", img->wd, img->ht);
  size_t off = 0;
  while((len + 1 + off) & 0xf) off++;
  while(off-- > 0) fprintf(f, "0");
  fprintf(f, "\n");
  fwrite(img->pixel, img->wd*img->ht*3, sizeof(float), f);
  fclose(f);
  return 0;
}

