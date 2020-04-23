#ifndef VMLT_MULTICHAIN_H
#define VMLT_MULTICHAIN_H

#include "points.h"
#include "pathspace/halfvec.h"
#include "pathspace/multichain.h"

void *multichain_init()
{
  return 0;
}

void multichain_cleanup(void *data) { }

float multichain_suitability(const path_t *p, void *data)
{
  // currently only implemented for paths started at the eye
  if(!(p->v[0].mode & s_sensor)) return 0.0f;
  return 1.0f;
}

float multichain_mutate(path_t *curr, path_t *tent, void *data)
{
  // TODO: optimise this copy away! apparently we need to set shader id somewhere:
  *tent = *curr;
  const int tid = common_get_threadid();
  double T = 1.0;
  float g1, g2;
  const float r1 = points_rand(rt.points, tid);
  const float r2 = points_rand(rt.points, tid);
  sample_gaussian(r1, r2, &g1, &g2);
  // const float stepsize = .4 * HALFVEC_MUTATION_STEP;
  // tent->sensor.pixel_i += stepsize * g1;
  // tent->sensor.pixel_j += stepsize * g2;
  tent->sensor.pixel_i += g1 * 5.f;
  tent->sensor.pixel_j += g2 * 5.f * view_width() / (float)view_height();
  // mutate point on aperture, after mutating outgoing direction (pixel)
  tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_set = 1;
  tent->sensor.pixel_set = 1;
  tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
  tent->time = curr->time; // <= TODO

  // camera_mutate_aperture(rt.cam, tent, points_rand(rt.points, tid), points_rand(rt.points, tid), 0.2f);
  // T /= camera_pdf_mutate_aperture(rt.cam, curr, tent, 0.2f);
  // T *= camera_pdf_mutate_aperture(rt.cam, tent, curr, 0.2f);
  T *= multichain_perturb(curr, tent, 0, curr->length-1);
  if(!(T > 0.0))
  {
    tent->length = 0;
    return 0.0f;
  }

  return T;
}

void multichain_print_info(FILE *f, void *data)
{
  fprintf(f, "         : multi-chain perturbation\n");
}

#endif
