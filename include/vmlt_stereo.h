#pragma once

#include "view.h"
#include "sampler_common.h"

void *stereo_init()
{
  return 0;
}

void stereo_cleanup(void *data) { }

float stereo_suitability(const path_t *p, void *data)
{
  // no need to reproject to the only camera:
  // if((p->v[0].mode & s_sensor) && (p->length > 2) && (view_pdf_camid(0) < 1.0f)) return 1.0f;
  // if((p->v[0].mode & s_sensor) && (p->length > 2))
  if(p->length > 1)
    return 1.0f; // also work on dof
  return 0.0f;
}

float stereo_mutate(path_t *curr, path_t *tent, void *data)
{
  // new independent sample, the render module already set the index
  const int tid = common_get_threadid();
  const int camid = view_sample_camid(points_rand(rt.points, tid));
  // if(camid == curr->sensor.camid) return 0.0f;
  *tent = *curr;

  tent->sensor.aperture_set = 1;
  tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.01f);
  tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.01f);
  tent->sensor.camid = camid;
  path_project(tent, 1, s_propagate_mutate);
  if(!primid_eq(tent->v[1].hit.prim, curr->v[1].hit.prim)) return 0.0f;
  for(int k=0;k<3;k++) // test for proximity, too:
    if(fabsf(tent->v[1].hit.x[k] - curr->v[1].hit.x[k]) >
        1e-4f*MAX(0.5, MAX(fabsf(tent->v[1].hit.x[k]), fabsf(curr->v[1].hit.x[k]))))
      return 0.0f;
  
  // evaluate only the first two segments
  const double f_curr = path_measurement_contribution_dx(curr, 0, MIN(curr->length-1, 2));
  const double f_tent = path_measurement_contribution_dx(tent, 0, MIN(curr->length-1, 2));
  // symmetric mutation, T_ct = T_tc
  return f_tent/f_curr;
}

void stereo_print_info(FILE *f, void *data)
{
  fprintf(f, "         : stereo mutation\n");
}
