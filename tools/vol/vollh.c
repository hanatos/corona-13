#include "vol/lighthierarchy.h"
#include "vol/trace.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
  srand48(666);
  if(argc < 2)
  {
    fprintf(stderr, "usage: vollh input.vol\n");
    exit(1);
  }
  vol_tree_t *tree = vol_open(argv[1]);
  if(!tree)
  {
    fprintf(stderr, "failed to open volume!\n");
    exit(1);
  }
#if 1
  const int lod = 0;
  const int N = 32;
  const float cullx[3] = {0, 0, 0};
  const float culln[3] = {1, 0, 0};
  for(int j=0;j<N;j++)
  for(int i=0;i<N;i++)
  for(int k=0;k<N;k++)
  {
    float pos[3] = {0.0f};
    float pdf = 0.0f;
    const float L = vol_lighthierarchy_sample_point(
        tree, 550.0f, 0.0f, lod, 0,
        drand48(), drand48(), drand48(),
        // (i+drand48())/(float)N,
        // (j+drand48())/(float)N,
        // (k+drand48())/(float)N,
        // 0, 0,
        cullx, culln,
        pos,
        &pdf);
    float pdf2 = vol_lighthierarchy_pdf_point(
        tree, 550.0f, 0.0f, lod, 0,
        // 0, 0,
        cullx, culln,
        pos);
    if (!(fabsf(pdf2 - pdf) < 1e-4f)) {
      fprintf(stderr, "pdf: %g pdf2: %g, pos: %g %g %g, L: %g\n", pdf, pdf2, pos[0], pos[1], pos[2], L);
      exit(0);
    }
    //assert(fabsf(pdf2 - pdf) < 1e-4f);
    fprintf(stdout, "%g %g %g %g\n", pos[0], pos[1], pos[2], L);
  }
#endif
  vol_close(tree);
  exit(0);
}
