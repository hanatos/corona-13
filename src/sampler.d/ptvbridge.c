#include "pathspace.h"
#include "pathspace/vbridge.h"
#include "pointsampler.h"

typedef struct sampler_t
{
  int keep;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) { }

static inline float
sampler_mis(path_t *path)
{
  double pdf2 = 0.0f;
  double pdf_pt = 1.0f;
  for(int k=0;k<path->length;k++)
  {
    pdf_pt *= path_pdf_extend(path, k);
    double pdf = 0.0f;
    double pdf_vb = 0.0f;
    if(k >= 1 && k < path->length-1)
      pdf_vb = vbridge_pdf(path, path->length-1, path->length-k-1);
    pdf = pdf_pt * pdf_vb;
    pdf2 += pdf*pdf;
  }
  pdf2 += pdf_pt * pdf_pt;
  double pdf = path_pdf(path);
  double w = (pdf*pdf)/pdf2;
  return w;
}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
      pointsampler_splat(path, path_throughput(path) * sampler_mis(path));
#if 0
      if(path->length > 3)
        if(path_russian_roulette(path, fminf(1.0, path->v[path->length-1].throughput/path->v[path->length-2].throughput)))
          return;
#endif
    }

    const int old_length = path->length;
    if(!vbridge_sample(path))
      pointsampler_splat(path, path_throughput(path) * sampler_mis(path));
    vbridge_pop(path, old_length);
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : pathtracer with verter bridges\n");
}
