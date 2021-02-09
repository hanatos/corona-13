#include "pathspace.h"
#include "pathspace/tech.h"
#include "pathspace/nee.h"
// XXX #include "pathspace/mvnee.h"
// #include "pathspace/mnee.h"

mf_t path_tech_vertex_pdf_as_sampled(path_t *p, int v)
{
  switch(p->v[v].tech)
  {
    case s_tech_extend:
      return path_pdf_extend(p, v);
    case s_tech_extend_adj:
      return path_pdf_extend_adjoint(p, v);
    case s_tech_nee:
      return nee_pdf(p, v);
    case s_tech_nee_adj:
      return nee_pdf_adjoint(p, v);
    case s_tech_mnee_chain:
      // by convention, we assemble all pdf in the end vertex of the chain (the
      // jacobian is global to the chain, individual pdfs make little sense).
      return mf_set1(1.0f);
    case s_tech_mnee_end:
      return mf_set1(0.0f);//mnee_pdf(p, v, 0); // XXX this is not doing what you expect it to.
    case s_tech_merged:
      assert(0);
      return mf_set1(0.0f/0.0f); // XXXphoton_pdf_path_merge(p, v);
      // XXX disabled for mf
    // case s_tech_mvnee:
      // return mvnee_pdf(p, v);
    default: // argh
      return mf_set1(0.0f);
  }
}

md_t path_tech_pdf_as_sampled(path_t *p)
{
  md_t pdf = md_set1(1.0);
  for(int k=0;k<p->length;k++)
    pdf = md_mul(pdf, mf_2d(path_tech_vertex_pdf_as_sampled(p, k)));
  return pdf;
}
