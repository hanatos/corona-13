#ifndef VMLT_VERTEXAREA_H
#define VMLT_VERTEXAREA_H

#include "prims.h"

typedef struct vmlt_vertexarea_t
{
  float *area;
  uint64_t num_prims;
}
vmlt_vertexarea_t;

void *vertexarea_init()
{
  vmlt_vertexarea_t *d = (vmlt_vertexarea_t *)malloc(sizeof(vmlt_vertexarea_t));
  d->num_prims = rt.prims->num_prims;
  d->area = (float *)malloc(sizeof(float)*d->num_prims);
  for(int k=0;k<rt.prims->num_prims;k++)
    d->area[k] = prims_get_area(rt.prims, rt.prims->primid[k]);
  for(int k=1;k<rt.prims->num_prims;k++)
    d->area[k] += d->area[k-1];
  for(int k=0;k<rt.prims->num_prims-1;k++)
    d->area[k] /= d->area[d->num_prims-1];
  d->area[d->num_prims-1] = 1.0f;
  return d;
}

void vertexarea_cleanup(void *data)
{
  vmlt_vertexarea_t *d = (vmlt_vertexarea_t *)data;
  free(d->area);
  free(d);
}

float vertexarea_suitability(const path_t *p, void *data)
{
  // can't construct path from scratch, so demand initialized sample:
  if(p->length) return 1.0f;
  return 0.0f;
}

float vertexarea_mutate(path_t *curr, path_t *tent, void *data)
{
  vmlt_vertexarea_t *d = (vmlt_vertexarea_t *)data;
  // take path curr, copy over, redistribute all vertices but the one on the eye and light
  memcpy(tent, curr, sizeof(path_t)); // wasteful, but whatever.
  float fullpdf = 1.0f; // account for volume mutations and such along the ray
  float pdf = 1.0f;
  const int tid = common_get_threadid();
  for(int v=1;v<tent->length-1;v++)
  {
    const float r0 = points_rand(rt.points, tid);
    const float r1 = points_rand(rt.points, tid);
    const float r2 = points_rand(rt.points, tid);
    uint64_t t = sample_cdf(d->area, d->num_prims, r0);
    // write to v[v].hit.x (and prim + uv, but actually we're only going to use hit.x)
    prims_sample_time(rt.prims, rt.prims->primid[t], r1, r2, &tent->v[v].hit, curr->time);
    if(path_project(tent, v, &pdf)) return 0.0f;
    fullpdf *= pdf;
  }
  // now add connect last vertex (use the old one)
  if(path_project(tent, tent->length-1, &pdf)) return 0.0f;
  fullpdf *= pdf;
  // could have done similar things by next event estimation:
  // tent->length--;
  // if(path_next_event(tent)) return 0.0f;

  // acceptance should be straight measurement contribution, mutation was already in vertex area measure and symmetric.
  const double f_curr = path_measurement_contribution_dx(curr, 0, curr->length-1);
  const double f_tent = path_measurement_contribution_dx(tent, 0, tent->length-1);
  // symmetric mutation, T_ct = T_tc
  return f_tent/f_curr;
}

void vertexarea_print_info(FILE *f, void *data)
{
  fprintf(f, "         : random vertex on geometry (debug mutation)\n");
}
#endif
