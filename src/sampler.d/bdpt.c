/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "sampler.h"
#include "pointsampler.h"
#include "pathspace/nee.h"

// mis bdpt

typedef struct configuration_t
{
  float sum; // total path contribution for all configurations
  float c[PATHSPACE_MAX_VERTS+1][PATHSPACE_MAX_VERTS+1]; // path contribution for [npt][nlt]
}
configuration_t;

static inline void _configuration_init(configuration_t *c)
{
  memset(c, 0, sizeof(configuration_t));
}

#if 0
static inline void _configuration_update(configuration_t *c)
{
  // update total sum of contributions:
  c->sum = 0.0f;
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
    for(int j=0;j<=PATHSPACE_MAX_VERTS;j++)
      c->sum += c->c[i][j];
}

static inline void _configuration_accum(configuration_t *c, int npt, int nlt, float contrib)
{
  if(!(contrib >= 0.0)) return;
  common_atomic_add(&c->c[npt][nlt], contrib);
}

static inline void _configuration_print(configuration_t *c)
{
  _configuration_update(c);
  fprintf(stderr, "# path configuration vs. contribution (eye down, light right):\n");
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
  {
    // fprintf(stderr, "[%f]  ", c->s[i]);
    for(int j=0;j<=PATHSPACE_MAX_VERTS;j++)
      fprintf(stderr, "%f  ", c->c[i][j]/c->sum);
    fprintf(stderr, "\n");
  }
}
#else
static inline void _configuration_accum(configuration_t *c, int npt, int nlt, mf_t contrib) {}
static inline void _configuration_print(configuration_t *c) {}
#endif

typedef struct sampler_t
{
  configuration_t conf;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  _configuration_init(&s->conf);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  _configuration_print(&s->conf);
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}

void sampler_clear(sampler_t *s)
{
  _configuration_init(&s->conf);
}

static inline void asseq(mf_t u, mf_t v)
{
  // assert(fabsf(u - v) < 5e-2f*MAX(fabsf(u), fabsf(v)));
}

// this returns the sum of the pdf of all techniques used
// to construct the given path. we want it in projected solid
// angle (dwp), which means we discard all occurrences of
// geometry terms. this is meaningful only for the full path,
// vertex probabilities will depend on tracing direction.
md_t sampler_sum_pdf_dwp(path_t *path)
{
  const int light_dir = (path->v[0].mode & s_emit);
  if(path->length == 2)
  { // only path tracing does this:
    if(light_dir)
      return mf_2d(path_pdf_extend_adjoint(path, 0) * path_pdf_extend_adjoint(path, 1) / path_G(path, 1));
    return mf_2d(path_pdf_extend(path, 0) * path_pdf_extend(path, 1) / path_G(path, 1));
  }

  md_t pdf = md_set1(0.0);

  // vertex area pdf as if shot from
  md_t vpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  md_t vpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  md_t vpdf_nee_fwd, vpdf_nee_adj; // next event estimation pdfs on both ends

  // evaluate pdfs and adjoint pdfs (adjoint is wrt to the current path, always the opposite direction)
  vpdf_nee_fwd = mf_2d(nee_pdf(path, path->length-1)/path_G(path, path->length-1));
  vpdf_nee_adj = mf_2d(nee_pdf_adjoint(path, 0)/path_G(path, 1));
  for(int k=0;k<path->length;k++)
  {
    vpdf_fwd[k] = k ? mf_2d(path_pdf_extend(path, k)/path_G(path, k)) : md_set1(1.0);
    vpdf_adj[k] = (k==path->length-1) ? md_set1(1.0) : mf_2d(path_pdf_extend_adjoint(path, k)/path_G(path, k+1));
  }
  for(int k=1;k<path->length;k++)    vpdf_fwd[k] = md_mul(vpdf_fwd[k], vpdf_fwd[k-1]);
  for(int k=path->length-2;k>=0;k--) vpdf_adj[k] = md_mul(vpdf_adj[k], vpdf_adj[k+1]);

  // use these cached pdfs to construct pdfs for all techniques
  // num_fwd is how many vertices were created from the same side as we were constructed
  for(int num_fwd=0;num_fwd<=path->length;num_fwd++)
  {
    // specular connections have zero probability (but pdf() will return non-zero)
    if((num_fwd > 0 && num_fwd < path->length) &&
       ((path->v[num_fwd-1].mode & s_specular) ||
        (path->v[num_fwd].mode & s_specular)))
      continue;

    md_t opdf = md_set1(0.0);
    if(num_fwd == 0)
    {
      // all vertices created from the other side.
      if(!light_dir) opdf = md_set1(0.0); // lt can't hit the lens
      else opdf = vpdf_adj[0]; // adj pt pdf
    }
    else if(num_fwd == 1)
    {
      // next event estmiation from the other side
      opdf = md_mul(vpdf_adj[1], vpdf_nee_adj);
      if(path->length < 3) opdf = md_set1(0.0); // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length-1)
    {
      // next event estimation from our side
      opdf = md_mul(vpdf_fwd[path->length-2], vpdf_nee_fwd);
      if(path->length < 3) opdf = md_set1(0.0); // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length)
    {
      // scatter and hit by chance. lt doesn't do that
      if(light_dir) opdf = md_set1(0.0);
      else opdf = vpdf_fwd[path->length-1]; // pt pdf
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      opdf = md_mul(vpdf_fwd[num_fwd-1], vpdf_adj[num_fwd]);
      
      opdf = md_div(opdf, md_set1(path_G(path, num_fwd)));
    }
    pdf = md_add(pdf, opdf);
  }

  return pdf;
}

// return mis weight
// light_v is the number of vertices created from the light (i.e. 0, 1, .. length-1 => pt ptdl .. lt)
static inline mf_t sampler_mis(path_t *path, int light_v)
{
  const int light_dir = (path->v[0].mode & s_emit);
  if(path->length == 2)
  {
    if(light_v == 0) return mf_set1(1.0f);
    else return mf_set1(0.0f); // no technique other than pt is good at sampling directly visible light sources
  }
  // got a path with l=light_v light vertices and e=path->length-light_v eye vertices.
  // if(!light_dir || light_v != path->length-1) return 0; // XXX lt
  // if(light_dir || light_v > 0) return 0; // XXX pt

  // get our pdf as sampled (still works with connected paths, is just a product of vertex area pdfs)
  md_t pdf = path_pdf(path);

  // vertex area pdf as if shot from
  mf_t vpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  mf_t vpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  mf_t vpdf_nee_fwd, vpdf_nee_adj; // next event estimation pdfs on both ends

  // evaluate pdfs and adjoint pdfs (adjoint is wrt to the current path, always the opposite direction)
  const int joint = light_dir ? light_v : (path->length-light_v);
  vpdf_nee_fwd = nee_pdf(path, path->length-1);
  vpdf_nee_adj = nee_pdf_adjoint(path, 0);
  for(int k=0;k<joint;k++)
  {
    vpdf_fwd[k] = path->v[k].pdf;
    vpdf_adj[k] = path_pdf_extend_adjoint(path, k);
    // FIXME: rough dielectric transmit for light tracing seems to have problems with this:
    asseq(vpdf_fwd[k], path_pdf_extend(path, k));
  }
  for(int k=joint;k<path->length;k++)
  {
    vpdf_fwd[k] = path_pdf_extend(path, k);
    if(joint < path->length-1)
    { // else v.pdf is nee pdf
      vpdf_adj[k] = path->v[k].pdf;
      asseq(vpdf_adj[k], path_pdf_extend_adjoint(path, k));
    }
    vpdf_adj[k] = path_pdf_extend_adjoint(path, k);
  }

  // TODO: port to mf_t
  // assert(vpdf_fwd[0] == 1.0);
  // assert(vpdf_adj[path->length-1] == 1.0);

  for(int k=1;k<path->length;k++)    vpdf_fwd[k] *= vpdf_fwd[k-1];
  for(int k=path->length-2;k>=0;k--) vpdf_adj[k] *= vpdf_adj[k+1];

  // use these cached pdfs to construct pdfs for all other techniques
  md_t sum_other_pdf2 = md_set1(0.0);
  // num_fwd is how many vertices were created from the same side as we were constructed
  for(int num_fwd=0;num_fwd<=path->length;num_fwd++)
  {
    // specular connections have zero probability (but pdf() will return non-zero)
    // only add our own pdf (num_fwd == joint) to improve numerical stability, see below.
    if((num_fwd != joint && num_fwd > 0 && num_fwd < path->length) &&
       ((path->v[num_fwd-1].mode & s_specular) ||
        (path->v[num_fwd].mode & s_specular)))
      continue;

    md_t opdf = md_set1(1.0);
    if(num_fwd == joint)
    {
      opdf = pdf; // that's us.
      // it's essential to sync these pdf, because if we compute it again here it might be off
      // by an order of magnitude because of floating point inaccuracies. we don't want to mess
      // with our overall weight normalisation..
#if 0
      // TODO: assert backward pdf is the same etc
      double pdf2 = vpdf_fwd[num_fwd-1];
      if(num_fwd < path->length-1)
        pdf2 *= vpdf_adj[num_fwd];
      else if(num_fwd == path->length-1)
        pdf2 *= vpdf_nee_fwd;
      asseq(opdf, pdf2);
      // TODO:
      // double pdf3 = 1.0;
      // path_t tmp;
      // path_reverse(&tmp, path);
      // ...
#endif
    }
    else if(num_fwd == 0)
    {
      // all vertices created from the other side.
      if(!light_dir) opdf = md_set1(0.0); // lt can't hit the lens
      else opdf = mf_2d(vpdf_adj[0]); // adj pt pdf
    }
    else if(num_fwd == 1)
    {
      // next event estmiation from the other side
      opdf = mf_2d(vpdf_adj[1] * vpdf_nee_adj);
      if(path->length < 3) opdf = md_set1(0.0); // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length-1)
    {
      // next event estimation from our side
      opdf = mf_2d(vpdf_fwd[path->length-2] * vpdf_nee_fwd);
      if(path->length < 3) opdf = md_set1(0.0); // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length)
    {
      // scatter and hit by chance. lt doesn't do that
      if(light_dir) opdf = md_set1(0.0);
      else opdf = mf_2d(vpdf_fwd[path->length-1]); // pt pdf
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      opdf = mf_2d(vpdf_fwd[num_fwd-1] * vpdf_adj[num_fwd]);
    }
    sum_other_pdf2 = md_add(sum_other_pdf2, md_mul(opdf, opdf));
  }

  return md_2f(md_div(md_mul(pdf, pdf), md_set1(md_hsum(sum_other_pdf2))));
}

void sampler_create_path(path_t *path)
{
  // lt path
  path->v[0].mode = s_emit; // start at the light
  path->v[0].rand_beg = s_dim_num_pt_beg + 1 + s_dim_num_nee; // leave first dimensions for path tracer
  while(1)
  {
    if(path_extend(path)) break;
    if(nee_sample(path)) break;
    const int v2 = path->length-1;
    if(mf_any(mf_gt(path->v[v2].throughput, mf_set1(0.0f))) && (path->v[v2].mode & s_sensor))
    {
      const mf_t weight = sampler_mis(path, path->length-1);
      _configuration_accum(&rt.sampler->conf, 1, path->length-1, mf_mul(path->v[v2].throughput, weight));
      pointsampler_splat(path, mf_mul(path->v[v2].throughput, weight));
    }
    path_pop(path);
    // interleave random dimensions between paths, leave space for pt:
    path->v[path->length-1].rand_cnt += s_dim_num_extend + s_dim_num_nee;
  }
  // XXX still okay even in case of break?
  // undo random number count adjustments:
  for(int k=1;k<path->length;k++)
    path->v[k].rand_cnt -= s_dim_num_extend + s_dim_num_nee;

  // eye path
  path_t eye_path;
  // same qmc sequence and same stereo camera
  path_init(&eye_path, path->index, path->sensor.camid);
  // instruct kelemen to mutate further from here:
  // eye_path.v[0].rand_beg = path->v[path->length-1].rand_beg + path->v[path->length-1].rand_cnt;
  eye_path.lambda = path->lambda;
  eye_path.time = path->time;
  while(1)
  {
    if(path_extend(&eye_path)) break;
    const int v = eye_path.length-1;
    if((eye_path.v[v].mode & s_emit) && mf_any(mf_gt(eye_path.v[v].shading.em, mf_set1(0.0f))))
    {
      const mf_t weight = sampler_mis(&eye_path, 0);
      _configuration_accum(&rt.sampler->conf, eye_path.length, 0, mf_mul(path_throughput(&eye_path), weight));
      pointsampler_splat(&eye_path, path_throughput(&eye_path) * weight);
    }

    if(nee_sample(&eye_path)) break;
    const int v2 = eye_path.length-1;
    if(mf_any(mf_gt(eye_path.v[v2].throughput, mf_set1(0.0f))) && (eye_path.v[v2].mode & s_emit))
    {
      const mf_t weight = sampler_mis(&eye_path, 1);
      _configuration_accum(&rt.sampler->conf, eye_path.length-1, 1, mf_mul(eye_path.v[v2].throughput, weight));
      pointsampler_splat(&eye_path, eye_path.v[v2].throughput * weight);
    }
    path_pop(&eye_path);

    // interleaved random numbers:
    const int v3 = eye_path.length-1;
    if(v3 > 1)
      eye_path.v[v3].rand_cnt += s_dim_num_extend + s_dim_num_nee;
    else
      eye_path.v[v3].rand_cnt += s_dim_num_lt_beg + 1 + s_dim_num_nee;
  }
  // undo random number count adjustments:
  eye_path.v[1].rand_cnt -= s_dim_num_lt_beg + 1 + s_dim_num_nee;
  for(int k=2;k<eye_path.length;k++)
    eye_path.v[k].rand_cnt -= s_dim_num_extend + s_dim_num_nee;

  // we now constructed paths of type #light vertices = 0, 1, and length-1.
  // now construct the rest using path_connect.

  path_t path2; // make a copy, we'll destroy it's contents. could work in-place, but not with klemen mlt
  memcpy(&path2, path, sizeof(path_t));

  // remember eye path random numbers and some meta-info for clever mutations.
  // this is quite brittle, kmlt mutate() assumes the eye path comes later (s_emit and v[0] trigger special cases)
  for(int k=0;k<eye_path.length && k < PATHSPACE_MAX_VERTS-path->length;k++)
    path->v[path->length+k] = eye_path.v[k];
  path->length = MIN(PATHSPACE_MAX_VERTS, path->length + eye_path.length);

  const int eye_path_length = eye_path.length;

  // l==0 would be using all eye path vertices, which is already handled (pt)
  // l==1 would be next event estimation, lt will always produce the two first vertices
  //      in one go, so 2 is the minimum:
  for(int l=path2.length;l>=2;l--)
  {
    // e==0 would be hitting the lens by chance, we don't have that technique.
    // e==1 would be light tracing, which is already considered in the connection above.
    for(int e=eye_path_length;e>=2;e--)
    {
      // this will destroy the contents of path2 after path2.v[l-1],
      // so we're slowly going near it from the back side
      path2.length = l;
      eye_path.length = e;
      const mf_t throughput = path_connect(&path2, &eye_path);
      if(mf_all(mf_lte(throughput, mf_set1(0.f)))) continue; // might fail due to excessive path length
      const mf_t weight = sampler_mis(&path2, l);
      _configuration_accum(&rt.sampler->conf, e, l, mf_mul(throughput, weight));
      pointsampler_splat(&path2, mf_mul(throughput, weight));
    }
  }
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : mis-weighted bi-directional path tracer\n");
}
