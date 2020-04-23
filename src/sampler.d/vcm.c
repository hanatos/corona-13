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

#include "threads.h"
#include "sampler.h"
#include "pointsampler.h"
#include "pathspace/nee.h"
#include "pathspace/photon.h"
#include "pathspace/manifold.h"

typedef struct sampler_t
{
  photon_map_t *pmap;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = calloc(1, sizeof(sampler_t));
  // s->pmap = photon_init(10*view_width()*view_height());
  s->pmap = photon_init(view_width()*view_height());
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  photon_cleanup(s->pmap);
  free(s);
}

static void *trace_photon_paths(void *arg)
{
  sampler_t *s = rt.sampler;
  path_t path;
  for(int k=0;k<s->pmap->max_paths/rt.num_threads;k++)
  {
    path_init(&path,
        k + (common_get_threadid()*s->pmap->max_paths)/rt.num_threads,
        0); // camera id doesn't matter here, leave at 0
    path.v[0].mode = s_emit; // start at the light
    while(1)
      if(path_extend(&path)) break;

    if(photon_record_path(
          rt.sampler->pmap,
          &path,
          path.v[path.length-1].throughput))
      return 0; // cache full. TODO: this may throw our throughput estimates off. should not happen
  }
  return 0;
}

void sampler_prepare_frame(sampler_t *s)
{
  // TODO: maybe do this only every once in a while:
  photon_clear(s->pmap);
  // shoot photons
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_task_init(rt.threads->task + k, &rt.threads->pool, trace_photon_paths, s);
  pthread_pool_wait(&rt.threads->pool);

  // build kd tree
  photon_build(s->pmap);
  // fprintf(stderr, "photon paths: %d/%d verts %d\n", s->pmap->num_paths, s->pmap->max_paths, s->pmap->num_verts);
}

void sampler_clear(sampler_t *s)
{
  photon_clear(s->pmap);
}

// TODO: remove v.pdf and replace by generic pdf cache?
// TODO: write generic version sampler_mis_weight for erpt!
// return mis weight
// light_v is the number of vertices created from the light (i.e. 0, 1, .. length-1 => pt ptdl .. lt)
// TODO: remove light_v argument and use vertex->tech interface instead!
static inline float sampler_mis(path_t *path, int light_v)
{
  const int light_dir = (path->v[0].mode & s_emit);
  if(path->length == 2)
  {
    if(light_v == 0) return 1.0f;
    else return 0.0f; // no technique other than pt is good at sampling directly visible light sources
  }

  // get our pdf as sampled (still works with connected paths, is just a product of vertex area pdfs)
  // TODO: use path_tech_pdf_as_sampled() and then it'll work with erpt!
  double pdf = path_pdf(path);
  // assert(pdf > 0.0); // happens every blue moon if the product drops to 0
  if(!(pdf > 0.0)) return 0.0;

  // vertex area pdf as if shot from
  float vpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  float vpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  float vpdf_nee_fwd, vpdf_nee_adj; // next event estimation pdfs on both ends

  // evaluate pdfs and adjoint pdfs (adjoint is wrt to the current path, always the opposite direction)
  const int joint = light_dir ? light_v : (path->length-light_v);
  vpdf_nee_fwd = nee_pdf(path, path->length-1);
  vpdf_nee_adj = nee_pdf_adjoint(path, 0);
  for(int k=0;k<path->length;k++)
  {
    vpdf_fwd[k] = path_pdf_extend(path, k);
    vpdf_adj[k] = path_pdf_extend_adjoint(path, k);
  }
  assert(vpdf_fwd[0] == 1.0);
  assert(vpdf_adj[path->length-1] == 1.0);
  for(int k=1;k<path->length;k++)    vpdf_fwd[k] *= vpdf_fwd[k-1];
  for(int k=path->length-2;k>=0;k--) vpdf_adj[k] *= vpdf_adj[k+1];

  // use these cached pdfs to construct pdfs for all other techniques
  double sum_other_pdf2 = 0.0f;
  // num_fwd is how many vertices were created from the same side as we were constructed
  for(int num_fwd=0;num_fwd<=path->length;num_fwd++)
  {
    // specular connections have zero probability (but pdf() will return non-zero)
    // only add our own pdf (num_fwd == joint) to improve numerical stability, see below.
    if((num_fwd != joint && num_fwd > 0 && num_fwd < path->length) &&
       ((path->v[num_fwd-1].mode & s_specular) ||
        (path->v[num_fwd].mode & s_specular)))
      continue;

    double opdf = 1.0f;
    if(num_fwd == 0)
    {
      // all vertices created from the other side.
      if(!light_dir) opdf = 0.0f; // lt can't hit the lens
      else opdf = vpdf_adj[0]; // adj pt pdf
    }
    else if(num_fwd == 1)
    {
      // next event estmiation from the other side
      opdf = vpdf_adj[1] * vpdf_nee_adj;
      if(path->length < 3) opdf = 0.0; // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length-1)
    {
      // next event estimation from our side
      opdf = vpdf_fwd[path->length-2] * vpdf_nee_fwd;
      if(path->length < 3) opdf = 0.0; // we don't do nee from the first vertex
    }
    else if(num_fwd == path->length)
    {
      // scatter and hit by chance. lt doesn't do that
      if(light_dir) opdf = 0.0f;
      else opdf = vpdf_fwd[path->length-1]; // pt pdf
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      opdf = vpdf_fwd[num_fwd-1] * vpdf_adj[num_fwd];
    }
    sum_other_pdf2 += opdf*opdf;
  }
  // add up vcm pdf for all potential merge verts:
  for(int k=1;k<path->length-1;k++)
  {
    // somewhat wasteful to recompute the _extend and _extend_adjoint pdfs here
    double merge_pdf = photon_pdf_path_merge(rt.sampler->pmap, path, k);
    merge_pdf *= vpdf_fwd[k-1] * vpdf_adj[k+1];
    sum_other_pdf2 += merge_pdf * merge_pdf;
  }

  // very sad to have to clip this rubbish.
  // unfortunately photon maps do weird things to paths.
  if(!(sum_other_pdf2 > 0.0)) return 0.0;
  double weight = pdf*pdf/sum_other_pdf2;
  if(!(weight > 0.0)) return 0.0;
  return CLAMP(weight, 0.0f, 1.0f);
}

void sampler_create_path(path_t *path)
{
  // lt path
  path->v[0].mode = s_emit; // start at the light
  path->v[0].rand_beg = s_dim_num_pt_beg + 1 + s_dim_num_nee; // leave first dimensions for path tracer
#if 1
  while(1)
  {
    if(path_extend(path)) break;
    if(nee_sample(path)) break;
    const int v2 = path->length-1;
    if(path->v[v2].throughput > 0.0f && (path->v[v2].mode & s_sensor))
    {
      const float weight = sampler_mis(path, path->length-1);
      pointsampler_splat(path, path->v[v2].throughput * weight);
    }
    path_pop(path);
    // interleave random dimensions between paths, leave space for pt:
    path->v[path->length-1].rand_cnt += s_dim_num_extend + s_dim_num_nee;
  }
  // XXX still okay even in case of break?
  // undo random number count adjustments:
  for(int k=1;k<path->length;k++)
    path->v[k].rand_cnt -= s_dim_num_extend + s_dim_num_nee;
#endif

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
#if 1
    if((eye_path.v[v].mode & s_emit) && (eye_path.v[v].shading.em > 0))
    {
      const float weight = sampler_mis(&eye_path, 0);
      pointsampler_splat(&eye_path, path_throughput(&eye_path) * weight);
      break; // light sources don't reflect
      // note that this does not work with multiple emitting vertices along
      // participating media such as fire. however, considering these
      // multiple contributions has problems in the context of MIS: we'd
      // need to marginalise the probability to pick up emission at these
      // vertices over all possible paths branching off this vertex.
    }
#endif
#if 1
    path_t tmp = eye_path;
    // need to initialise vertex derivatives for photon lookup
    manifold_compute_derivatives(&eye_path, v);
    // TODO: need to make sure these are inited on the whole path for mis!
    if(!photon_path_merge(rt.sampler->pmap, &tmp))
    {
      // we could claim one more
      const float weight = sampler_mis(&tmp, tmp.length-eye_path.length);
      pointsampler_splat(&tmp, tmp.throughput * weight);
    }
#endif
#if 1
    if(nee_sample(&eye_path)) break;
    const int v2 = eye_path.length-1;
    if(eye_path.v[v2].throughput > 0.0f && (eye_path.v[v2].mode & s_emit))
    {
      const float weight = sampler_mis(&eye_path, 1);
      pointsampler_splat(&eye_path, eye_path.v[v2].throughput * weight);
    }
    path_pop(&eye_path);
#endif

    // interleaved random numbers:
    const int v3 = eye_path.length-1;
    if(v3 > 1)
      eye_path.v[v3].rand_cnt += s_dim_num_extend + s_dim_num_nee;
    else
      eye_path.v[v3].rand_cnt += s_dim_num_lt_beg + 1 + s_dim_num_nee;
  }
#if 1
  // undo random number count adjustments:
  eye_path.v[1].rand_cnt -= s_dim_num_lt_beg + 1 + s_dim_num_nee;
  for(int k=2;k<eye_path.length;k++)
    eye_path.v[k].rand_cnt -= s_dim_num_extend + s_dim_num_nee;

  // we now constructed paths of type #light vertices = 0, 1, and length-1.
  // now construct the rest using path_connect.

  path_t path2; // make a copy, we'll destroy its contents. could work in-place, but not with klemen mlt
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
      const float throughput = path_connect(&path2, &eye_path);
      if(throughput <= 0.f) continue; // might fail due to excessive path length
      const float weight = sampler_mis(&path2, l);
      pointsampler_splat(&path2, throughput * weight);
    }
  }
#endif
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : vertex connection and merging\n");
}
