#include "pfm.h"
#include "fft_solver.h"
#include "cgsolve.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main (int argc, char *argv[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <basename>\n", argv[0]);
    exit(1);
  }
  char fname[256];
  pfm_t img, grad_x, grad_y;
  int ret = 0;
  snprintf(fname, sizeof(fname), "%s.pfm", argv[1]);
  ret += pfm_load(&img, fname);
  snprintf(fname, sizeof(fname), "%s_grad_x.pfm", argv[1]);
  ret += pfm_load(&grad_x, fname);
  snprintf(fname, sizeof(fname), "%s_grad_y.pfm", argv[1]);
  ret += pfm_load(&grad_y, fname);
  if(ret)
  {
    fprintf(stderr, "could not load %s.pfm %s_grad_x.pfm %s_grad_y.pfm\n", argv[1], argv[1], argv[1]);
    exit(1);
  }

  // jaakko uses alpha=0.2 in his paper, but the plots seem to indicate 0.3..0.4 are better for more complicated scenes:
  // fourier_solve(img.pixel, grad_x.pixel, grad_y.pixel, img.wd, img.ht, .2f);
  cg_solve(img.pixel, grad_x.pixel, grad_y.pixel, img.wd, img.ht, .2f);

  pfm_write(&img, "reconstructed.pfm");
  exit(0);
}
