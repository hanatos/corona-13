#ifndef SAMPLER_PT_MNEE_H
#define SAMPLER_PT_MNEE_H

#include "pathspace.h"
#include "pathspace/mnee.h"
#include "pointsampler.h"

// std backward pathtracer with stubborn next event estimation

typedef struct sampler_t
{
  halfvec_stats_t *stats;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  s->stats = (halfvec_stats_t *)malloc(sizeof(halfvec_stats_t)*rt.num_threads);
  memset(s->stats, 0, sizeof(halfvec_stats_t)*rt.num_threads);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  free(s->stats);
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s)
{
  memset(s->stats, 0, sizeof(halfvec_stats_t)*rt.num_threads);
}

static inline float sampler_mis(path_t *path)
{
  halfvec_stats_t *stats = rt.sampler->stats + common_get_threadid();
#ifdef PTMNEE_BIAS // bias: kill paths that we can't create with mnee for stupid reasons (but should actually be able to)
  int transmit = 0;
  // mnee is possible if there exists a transmissive postfix in the path..
  // count length of maximum postfix of transmissive inner vertices
  for(int k=path->length-2;k>=0;k--)
    if(path->v[k].mode & s_transmit) transmit++;
    else break;
  // // .. that is connectable through a non-singular bsdf
  // if(transmit && !(path->v[path->length-2-transmit].material_modes & (s_diffuse | s_glossy)))
  //   transmit = 0;
#endif
  // compute sum of pdfs of all possible ways to create this path.
  // we trust pdf to be one of the pdf we recompute here for simplicity.
  double pdf2 = 0.0f;
  // 1) std backward pt, scattering only
  // 2) all possible postfixes as created by mnee
  double pdf_pt = 1.0f;

#ifdef PTMNEE_BIAS
  // count total mnee pdf in cases where it should have worked:
  double sum_pdf_mnee = 0.0;
#endif
  for(int k=0;k<path->length;k++)
  {
    pdf_pt *= path_pdf_extend(path, k);
    double pdf = 0.0f;
    // returns pdf of creating the postfix [k+1..end] via mnee
    double pdf_mnee = 0.0f;
    if(k >= 1 && k < path->length-1)
      pdf_mnee = mnee_pdf(path, k+1, stats);
#ifdef PTMNEE_BIAS
    // if connecting vertex is non-singular and we have a purely transmissive postfix
    if((k >= path->length-2-transmit) && // purely transmissive postfix from here to the end
       (k < path->length-1) && // not the last vertex, at least one left for nee
       (path->v[k].material_modes & (s_diffuse | s_glossy))) // connectable
      sum_pdf_mnee += pdf_mnee;
#endif
    pdf = pdf_pt * pdf_mnee;
    pdf2 += pdf*pdf;
  }
#ifdef PTMNEE_BIAS
  // bias: kill unlucky paths where mnee failed due to unknown reasons:
  if(sum_pdf_mnee == 0.0 && transmit) return 0.0;
  // if(pdf2 == 0.0 && transmit) return 0.0;
#else
  pdf2 += pdf_pt * pdf_pt;
#endif
  double pdf = path_pdf(path);
  double w = (pdf*pdf)/pdf2;
  // XXX if(!(w > 0.0f) || !(w < FLT_MAX)) return 0.0f;
  return MIN(1.0f, w);
}

void sampler_create_path(path_t *path)
{
  halfvec_stats_t *stats = rt.sampler->stats + common_get_threadid();
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
#ifdef PTMNEE_BIAS // pretend we tagged the surface, so we don't need any mis.
      // if(path->length <= 2)
#endif
      pointsampler_splat(path, path_throughput(path) * sampler_mis(path));
      if(path->length > 3)
        if(path_russian_roulette(path, fminf(1.0, path->v[path->length-1].throughput/path->v[path->length-2].throughput)))
          return;
    }

    const int old_length = path->length;
    if(!mnee_sample(path, stats))
      pointsampler_splat(path, path_throughput(path) * sampler_mis(path));
    mnee_pop(path, old_length);
  }
}

static inline void mnee_stats_accum(halfvec_stats_t *accum, const halfvec_stats_t *stats)
{
  memset(accum, 0, sizeof(halfvec_stats_t));
  for(int k=0;k<rt.num_threads;k++)
  {
    accum->newton_walks += stats[k].newton_walks;
    accum->successful_newton_walks += stats[k].successful_newton_walks;
    accum->project_failed += stats[k].project_failed;
    accum->propagate_failed += stats[k].propagate_failed;
    accum->perturb_failed += stats[k].perturb_failed;
    accum->no_throughput += stats[k].no_throughput;
    accum->mutations += stats[k].mutations;
    accum->mutations_proposed += stats[k].mutations_proposed;
    accum->reverse_check_failed += stats[k].reverse_check_failed;
    accum->sum_iterations += stats[k].sum_iterations;
    accum->sum_successful_reduces += stats[k].sum_successful_reduces;
    accum->raydiff_compute_failed += stats[k].raydiff_compute_failed;
    accum->error_increased += stats[k].error_increased;
  }
}

void sampler_print_info(FILE *fd)
{
  halfvec_stats_t stats;
  mnee_stats_accum(&stats, rt.sampler->stats);
  fprintf(fd, "sampler  : pathtracer with stubborn next event estimation and mis\n");
  fprintf(fd, "           successful walks       %.02f%%\n", 100.0f*stats.successful_newton_walks/(float)stats.newton_walks);
  fprintf(fd, "           average iterations     %.02f\n", stats.sum_iterations/(float)stats.successful_newton_walks);
  fprintf(fd, "           successful step reduce %.02f%%\n", 100.0f*stats.sum_successful_reduces/(float)stats.sum_iterations);
  fprintf(fd, "           error increased        %.02f%%\n", 100.0f*stats.error_increased/(float)stats.sum_iterations);
}

#endif
