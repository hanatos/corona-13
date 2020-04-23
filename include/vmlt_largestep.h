#pragma once

#include "pathspace/vmlt.h"
#include "sampler.h"

void *largestep_init()
{
  return 0;
}

void largestep_cleanup(void *data) { }

float largestep_suitability(const path_t *p, void *data)
{
  return 1.0f;
}

float largestep_mutate(path_t *curr, path_t *tent, void *data)
{
  // new independent sample, the render module already set the index
  path_init(tent, tent->index, tent->sensor.camid);

  sampler_create_path(tent);
  const double f_curr = sampler_throughput(curr);
  // tent could use the faster path_throughput() function instead, but this will not consider all possible
  // ways to construct the path, it will be different for bdpt1 (f/sum(pdf) vs f/pdf)
  // and matches spot-on for pt:
  // if(path_throughput(tent) > 0.0)
  //   fprintf(stderr, "tent throughput %g %g %g\n", path_throughput(tent), sampler_throughput(tent), path_measurement_contribution_dx(tent)/path_pdf(tent));
  const double f_tent = sampler_throughput(tent);
  if(f_curr <= 0.0 && f_tent > 0.0) return 1.0;
  // symmetric mutation, T_ct = T_tc
  return f_tent/f_curr;
}

void largestep_print_info(FILE *f, void *data)
{
  fprintf(f, "         : large step mutation\n");
}
