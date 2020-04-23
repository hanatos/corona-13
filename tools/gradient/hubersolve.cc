#include "hubersolve.hh"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct img_t
{
  float *pixel;
  int64_t wd, ht;
}
img_t; 

int load_pfm(img_t *img, const char *filename)
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

int write_pfm(const img_t *img, const char *filename)
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


int main (int argc, char *argv[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <basename>\n", argv[0]);
    return 1;
  }
  char fname[256];
  img_t img, grad_x, grad_y;
  int ret = 0;
  snprintf(fname, sizeof(fname), "%s.pfm", argv[1]);
  ret += load_pfm(&img, fname);
  snprintf(fname, sizeof(fname), "%s_grad_x.pfm", argv[1]);
  ret += load_pfm(&grad_x, fname);
  snprintf(fname, sizeof(fname), "%s_grad_y.pfm", argv[1]);
  ret += load_pfm(&grad_y, fname);
  if(ret)
  {
    fprintf(stderr, "could not load %s.pfm %s_grad_x.pfm %s_grad_y.pfm\n", argv[1], argv[1], argv[1]);
    return 1;
  }

  // jaakko uses alpha=0.2 in his paper, but the plots seem to indicate 0.3..0.4 are better for more complicated scenes:
  huber_solve(img.pixel, grad_x.pixel, grad_y.pixel, img.wd, img.ht, .1f);

  write_pfm(&img, "reconstructed.pfm");
  return 0;
}
