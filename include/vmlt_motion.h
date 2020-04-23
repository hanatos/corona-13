#pragma once

#include "vmlt.h"
#include "view.h"
#include "prims.h"
#include <float.h>

void *motion_init()
{
  return 0;
}

void motion_cleanup(void *data) { }

float motion_suitability(const path_t *p, void *data)
{
  if((p->v[0].mode & s_sensor) && (p->length > 2)) return 1.0f;
  // TODO: check if path has primitives with motion blur at all
  return 0.0f;
}

float motion_mutate(path_t *curr, path_t *tent, void *data)
{
  // new independent sample, the render module already set the index
  *tent = *curr;

  // sample time (for some reason this makes the image go black after a while, needs debugging)
  // TODO: do small mutations instead of completely independent random?
  tent->time = view_sample_time(tent, points_rand(rt.points, common_get_threadid()));
  for(int v=1;v<tent->length;v++)
  {
    if(!primid_invalid(tent->v[v].hit.prim)) // TODO: check what this does for volumes.
      prims_retime(rt.prims, tent->v[v].hit.prim, &tent->v[v].hit, tent->time);
    float x[3] = {tent->v[v].hit.x[0], tent->v[v].hit.x[1], tent->v[v].hit.x[2]};
    tent->e[v].dist = FLT_MAX; // TODO: keep that distance for volumes?
    if(path_project(tent, v, s_propagate_mutate)) return 0.0f;
    if(!primid_eq(tent->v[v].hit.prim, curr->v[v].hit.prim)) return 0.0f;
    for(int k=0;k<3;k++) // test for proximity, too:
      if(tent->v[v].hit.x[k] - x[k] >
          1e-4f*MAX(fabsf(tent->v[v].hit.x[k]), fabsf(x[k])))
        return 0.0f;
  }
  
  const double f_curr = path_measurement_contribution_dx(curr, 0, curr->length-1);
  const double f_tent = path_measurement_contribution_dx(tent, 0, tent->length-1);
  // FIXME: i think this is a lie (needs to be transformed for G terms)
  // symmetric mutation, T_ct = T_tc
  return f_tent/f_curr;
}

void motion_print_info(FILE *f, void *data)
{
  fprintf(f, "         : motion blur mutation\n");
}
