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

static inline mf_t
sampler_mis(path_t *path)
{
  md_t our = path_pdf(path);
  md_t sum = md_set1(0.0);
  md_t pdf_prefix = md_set1(1.0);
  for(int k=0;k<path->length;k++)
  {
    pdf_prefix = md_mul(pdf_prefix, mf_2d(path_pdf_extend(path, k)));
    md_t pdf_vb = md_set1(0.0);
    if(k >= 1 && k < path->length-1)
      pdf_vb = mf_2d(vbridge_pdf(path, path->length-1, path->length-k-1));
    md_t pdf = md_mul(pdf_prefix, pdf_vb);
    sum = md_add(sum, pdf);
    if(k == path->length-1)
      sum = md_add(sum, pdf_prefix); // plain pt technique
  }
  return mf_div(md_2f(our), mf_set1(mf_hsum(md_2f(sum))));
}

void sampler_create_path(path_t *path)
{
  while(1)
  {
    if(path_extend(path)) return;
    const int v = path->length-1;
    if(path->v[v].mode & s_emit)
    {
      // XXX DEBUG pointsampler_splat(path, path_throughput(path) * sampler_mis(path));
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
  fprintf(fd, "sampler  : pathtracer with vertex bridges\n");
}
