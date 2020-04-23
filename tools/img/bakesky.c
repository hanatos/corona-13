#include "../../src/shaders/daylight.h"
#include "screenshot.h"
#include "corona_common.h"
#include "filter.h"
#include "pathspace.h"
#include "spectrum.h"

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  const int wd = 4096;
  const int ht = 2048;
  float *buf = (float *)malloc(sizeof(float)*3*wd*ht);
  memset(buf, 0, 3*wd*ht*sizeof(float));
  void *data;
  fprintf(stderr, "please enter shader def line including new line, -z up (ex: 0 -1 -1 3 # \\n ^d)\n");
  sky_daylight_init(stdin, &data);
  path_t path;
  memset(&path, 0, sizeof(path_t));
  path.length = 2;
  path.time = 0.;
  int progress = 0;
#pragma omp parallel for schedule(static) default(shared) firstprivate(path)
  for(int j=0;j<ht;j++)
  {
    if((progress & 7) == 7)
      fprintf(stderr, "line %d/%d\r", progress, ht);
#pragma omp atomic
    progress++;
    for(int i=0;i<wd;i++)
    {
      const float phi = 2.0*M_PI*i/(float)wd;
      const float theta = M_PI*j/(float)ht;
      path.e[1].omega[0] = cosf(phi) * sinf(theta);
      path.e[1].omega[1] = sinf(phi) * sinf(theta);
      path.e[1].omega[2] = cosf(theta);
      const int numl = 20;
      for(int l=0;l<numl;l++)
      {
        path.lambda = 400.0 + l/(numl-1.0)*300.0;
        float xyz[3];
        spectrum_p_to_xyz(path.lambda, sky_daylight(&path, data), xyz);
        for(int k=0;k<3;k++) xyz[k] /= (float)numl;

        filter_splat(i, j, xyz, buf, 3, 0, 3, wd, ht);
        // for(int k=0;k<3;k++) buf[3*(i + wd*j) + k] += xyz[k];
      }
    }
  }
  fprintf(stderr, "done.                     \n");
  screenshot_write("sky", buf, 3, wd, ht);
  free(buf);
  free(data);
  exit(0);
}
