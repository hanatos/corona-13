#pragma once

// keep track of sampling technique which was used to construct a particular vertex:

typedef enum path_tech_t
{
  s_tech_none       = 0,  // 0 uninited
  s_tech_extend     = 1,  // extended via regular pt
  s_tech_extend_adj = 2,  // extended from other side
  s_tech_nee        = 3,  // next event estimation
  s_tech_nee_adj    = 4,  // next event estimation from other side
  s_tech_mnee_chain = 5,  // inner vertex in specular chain of mnee
  s_tech_mnee_end   = 6,  // end vertex of specular chain of mnee
  s_tech_merged     = 7,  // sampled from both ends and merged via photon map
  s_tech_mvnee      = 8,  // multiple vertex next event estimation
  s_tech_equiangular= 9,  // equiangular sampling for nee in volumes
}
path_tech_t;

// return the pdf of this vertex as it has been sampled, in vertex area
// measure.  (this should mostly be equivalent to the cached v.pdf, but not
// always: for instance if the path is used in an MLT mutation and we're
// interested in the MC equivalent).
mf_t path_tech_vertex_pdf_as_sampled(path_t *p, int v);

// same but for the full path
md_t path_tech_pdf_as_sampled(path_t *p);
