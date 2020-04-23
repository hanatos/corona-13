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
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "sampler.h"
#include "pointsampler.h"

// bdpt with just one connection. useful for metropolis.
// the number of eye path and light path vertices needs to be sampled from the outside and
// stored on the path.

//=======================================================================================================
typedef struct configuration_t
{
  float    sum; // total path contribution for all configurations
  float    c[PATHSPACE_MAX_VERTS+1][PATHSPACE_MAX_VERTS+1]; // path contribution for [npt][nlt]
  float    s[PATHSPACE_MAX_VERTS+1];                        // boundary sums [npt] (nlt summed away), normalised to sum to 1
  uint64_t n[PATHSPACE_MAX_VERTS+1][PATHSPACE_MAX_VERTS+1]; // number of samples to estimate c[][] and p[][]
}
configuration_t;

static inline void _configuration_update(configuration_t *c)
{
  // update boundary sums and total sum of contributions:
  c->sum = 0.0f;
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
  {
    c->s[i] = 0.0f;
    for(int j=0;j<=PATHSPACE_MAX_VERTS;j++)
      c->s[i] += c->c[i][j];
    c->sum += c->s[i];
  }
  // start at 1, 0 eye vertices doesn't exist (no by-chance lens intersection)
  for(int i=1;i<=PATHSPACE_MAX_VERTS;i++)
    c->s[i] /= c->sum;
}

static inline float _configuration_p(configuration_t *cc, int npt, int nlt)
{
  configuration_t *c = cc+common_get_threadid();
  return c->c[npt][nlt]/c->sum;
}

static inline void _configuration_accum(configuration_t *cc, int npt, int nlt, float contrib)
{
  // TODO: need to call configuration_accum in case of 0 contributions!
  if(!(contrib >= 0.0)) return;
  configuration_t *c = cc+common_get_threadid();
  c->c[npt][nlt] = c->c[npt][nlt]*((float)c->n[npt][nlt]/(c->n[npt][nlt]+1.0f)) + contrib/(c->n[npt][nlt]+1.0f);
  c->n[npt][nlt]++;
  // TODO: do selective update only on npt-th row!
  _configuration_update(c);
}

static inline void _configuration_print(const configuration_t *c)
{
  fprintf(stderr, "# path configuration vs. contribution (eye down, light right):\n");
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
  {
    // fprintf(stderr, "[%f]  ", c->s[i]);
    for(int j=0;j<=PATHSPACE_MAX_VERTS;j++)
      fprintf(stderr, "%f  ", c->c[i][j]/c->sum);
    fprintf(stderr, "\n");
  }
}

static inline void _configuration_init(configuration_t *c)
{
  // init with some empric, half way sane defaults:
  // set number of samples n == 100, so these will be overwritten smoothly with real samples,
  // with a half life of 100 samples.
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
  { // for all #pt verts
    for(int j=0;j<=PATHSPACE_MAX_VERTS;j++)
    { // for all #lt verts
      const float b = 100.0f; // TODO needs to depend on mean image brightness! which we don't know! so be conservative and use high b
      c->c[i][j] = b * fmaxf(1e-3f, expf(-(i+j)*(i+j)/20.0f));
      if(j == 1 || i == 1 || j == 0)
        c->c[i][j] = fmaxf(c->c[i][j], b*.1f); // pt and next event estimation is usually good
      c->n[i][j] = 1000;
      if(i+j > PATHSPACE_MAX_VERTS) c->c[i][j] = 0.0f;
    }
  }
  // kill meaningless combinations:
  for(int i=0;i<=PATHSPACE_MAX_VERTS;i++)
    c->c[0][i] = 0.0f; // no technique lt to intersect lens by chance
  c->c[1][0] = 0.0f; // need at least two verts
  c->c[1][1] = 0.0f; // can't both have only one vertex (means next event estimation)
  _configuration_update(c);
}

// TODO: helper function somewhere global? same with mis weights?
static inline void _configuration_sample(path_t *path, const configuration_t *cc)
{
  const int tid = common_get_threadid();
  const configuration_t *c = cc + tid;
  const float r1 = points_rand(rt.points, tid);
  const float r2 = points_rand(rt.points, tid);
  // sample 2d distribution for pt and lt vertex counts, do this via linear search.
  // building a cdf would be a linear scan, too, and those p change often (and aren't very large)
  int npt = 0;
  for(float cdf = c->s[0];
      npt < PATHSPACE_MAX_VERTS && (cdf < r1 || c->s[npt] == 0);
      npt++)
    cdf += c->s[npt+1];
  int nlt = 0;
  for(float cdf = c->c[npt][nlt]/(c->sum*c->s[npt]);
      nlt < PATHSPACE_MAX_VERTS-npt && (cdf < r2 || c->c[npt][nlt] == 0.0);
      nlt++)
    cdf += c->c[npt][nlt+1]/(c->sum*c->s[npt]);

  path->mmlt.pt_verts = npt;
  path->mmlt.lt_verts = nlt;
  path->mmlt.p = c->c[npt][nlt]/c->sum;

  // _configuration_print(c);
  // fprintf(stderr, "path length: %d %d = %d, %g\n", npt, nlt, npt+nlt, path->mmlt.p);

  assert(path->mmlt.p > 0.0);
  assert(npt + nlt <= PATHSPACE_MAX_VERTS);
  assert(npt + nlt >= 2);
  assert(npt >= 1); // no such technique as npt == 0 (intersect camera pupil by chance from the light)
  assert(nlt >= 0);
}
//=======================================================================================================

typedef struct sampler_t
{
  configuration_t *conf;
}
sampler_t;

sampler_t *sampler_init()
{
  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  s->conf = (configuration_t *)malloc(sizeof(configuration_t)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
    _configuration_init(s->conf + k);
  return s;
}

void sampler_cleanup(sampler_t *s)
{
  // for(int k=0;k<rt.num_threads;k++)
    // _configuration_print(s->conf+k);
  free(s->conf);
  free(s);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s)
{
  for(int k=0;k<rt.num_threads;k++)
    _configuration_init(s->conf + k);
}

// return mis weight (same as bdpt proper with all connections)
// light_v is the number of vertices created from the light (i.e. 0, 1, .. length-1 => pt ptdl .. lt)
static inline float sampler_mis(path_t *path, int light_v)
{
  const int light_dir = (path->v[0].mode & s_emit);
  if(path->length == 2)
  {
    if(light_v == 0) return 1.0f;
    else return 0.0f; // no technique other than pt is good at sampling directly visible light sources
  }
  // got a path with l=light_v light vertices and e=path->length-light_v eye vertices.

  // get our pdf as sampled (still works with connected paths, is just a product of vertex area pdfs)
  double pdf = path_pdf(path);

  // vertex area pdf as if shot from
  float vpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  float vpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  float vpdf_nee_fwd, vpdf_nee_adj; // next event estimation pdfs on both ends

  // evaluate pdfs and adjoint pdfs (adjoint is wrt to the current path, always the opposite direction)
  const int joint = light_dir ? light_v : (path->length-light_v);
  vpdf_nee_fwd = path_pdf_next_event(path, path->length-1);
  vpdf_nee_adj = path_pdf_next_event_adjoint(path, 0);
  for(int k=0;k<joint;k++)
  {
    vpdf_fwd[k] = path->v[k].pdf;
    vpdf_adj[k] = path_pdf_extend_adjoint(path, k);
  }
  for(int k=joint;k<path->length;k++)
  {
    vpdf_fwd[k] = path_pdf_extend(path, k);
    vpdf_adj[k] = path->v[k].pdf;
  }

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
    if(num_fwd == joint)
    {
      opdf = pdf; // that's us.
      // it's essential to sync these pdf, because if we compute it again here it might be off
      // by an order of magnitude because of floating point inaccuracies. we don't want to mess
      // with our overall weight normalisation..
    }
    else if(num_fwd == 0)
    {
      // all vertices created from the other side.
      if(!light_dir) opdf = 0.0f; // lt can't hit the lens
      else for(int k=0;k<path->length-1;k++) opdf *= vpdf_adj[k]; // pt, last two created in one path_extend() call
    }
    else if(num_fwd == 1)
    {
      // next event estmiation from the other side
      for(int k=1;k<path->length-1;k++) opdf *= vpdf_adj[k]; // len-2 and len-1 created in one call again
      opdf *= vpdf_nee_adj;
    }
    else if(num_fwd == path->length-1)
    {
      // next event estimation from our side
      for(int k=1;k<num_fwd;k++) opdf *= vpdf_fwd[k]; // start at 1, first two coupled.
      opdf *= vpdf_nee_fwd;
    }
    else if(num_fwd == path->length)
    {
      // scatter and hit by chance. lt doesn't do that
      if(light_dir) opdf = 0.0f;
      else for(int k=1;k<path->length;k++) opdf *= vpdf_fwd[k]; // pt, first two created in one call
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      for(int k=1;k<num_fwd;k++) opdf *= vpdf_fwd[k];
      for(int k=num_fwd;k<path->length-1;k++) opdf *= vpdf_adj[k];
    }
    sum_other_pdf2 += opdf*opdf;
  }

  return pdf*pdf/sum_other_pdf2;
}


void sampler_create_path(path_t *path)
{
  path_t lt_path;
  path_init(&lt_path);
  lt_path.index = path->index;
  // fake random number beg to include first dims for pt/ptdl/lt sampling, so kmlt can pick it up.
  // XXX needs more love (=> mmlt)
  // lt_path.v[0].rand_beg = 2;

  // if nobody configured us from the outside, do it now:
  if(path->mmlt.pt_verts + path->mmlt.lt_verts == 0)
    _configuration_sample(path, rt.sampler->conf);

  // read mmlt style path->mmlt.pt_verts and path->mmlt.lt_verts for configuration
  const int npt = path->mmlt.pt_verts;
  const int nlt = path->mmlt.lt_verts;
  assert(npt >= 1); // no such technique as npt == 0 (intersect camera pupil by chance from the light)
  assert(nlt >= 0);
  assert(npt + nlt >= 2);
  assert(npt + nlt <= PATHSPACE_MAX_VERTS);
  assert(path->mmlt.p > 0);

  float weight = 0.0f;

  // lt path
  lt_path.v[0].mode = s_emit; // start at the light
  if(nlt > 1) // only 1 would be pt + next event estimation, handled below.
  while(lt_path.length < nlt)
  {
    // try to extend the light path. if there is no more energy
    // transported by it or it escapes to the envmap, we reject the sample.
    if(path_extend(&lt_path)) return;
    const int v = lt_path.length-1;
    if(lt_path.v[v].throughput <= 0.0f) return;
    if(lt_path.v[v].flags & s_environment) return;

    if(lt_path.length == nlt && npt == 1)
    { // pure light tracing case
      if(path_next_event(&lt_path)) return;
      const int v2 = lt_path.length-1;
      if(lt_path.v[v2].throughput > 0.0f && (lt_path.v[v2].mode & s_sensor))
      {
        weight = sampler_mis(&lt_path, lt_path.length-1);
        pointsampler_splat(&lt_path, lt_path.v[v2].throughput * weight/path->mmlt.p);
        _configuration_accum(rt.sampler->conf, npt, nlt, lt_path.v[v2].throughput * weight/path->mmlt.p);
        // let mlt know about our cool path:
        path_reverse(path, &lt_path);
        // at this point we're done, there is no eye path (npt == 1 is next event estimation)
        goto success;
      }
      return;
    }
  }
  assert(nlt == 1 || lt_path.length == nlt); // or else we should have rejected the sample
  assert(npt > 1); // or else we should really have handled that above.

  // eye path
  // instruct kelemen to mutate further from here:
  if(lt_path.length)
    path->v[0].rand_beg = lt_path.v[lt_path.length-1].rand_beg + lt_path.v[lt_path.length-1].rand_cnt;
  else path->v[0].rand_beg = 0;
  path->lambda = lt_path.lambda;
  path->time = lt_path.time;
  while(path->length < npt)
  {
    if(path_extend(path)) goto no_path;
    const int v = path->length-1;
    if(path->v[v].throughput <= 0.0f) goto no_path;
    // pt case
    // TODO: should we discard those, too? seems very wasteful
    if(path->length == npt && nlt == 0)
    if(path->v[v].mode & s_emit)
    {
      weight = sampler_mis(path, 0);
      const float throughput = path_throughput(path);
      if(throughput > 0.0)
      {
        pointsampler_splat(path, throughput * weight / path->mmlt.p);
        _configuration_accum(rt.sampler->conf, npt, nlt, throughput * weight/path->mmlt.p);
        goto success;
      }
      goto no_path;
    }
    // envmap emission will be treated above already
    if(path->v[v].flags & s_environment) goto no_path;

    if(path->length == npt && nlt == 1)
    {
      if(path_next_event(path)) goto no_path;
      const int v2 = path->length-1;
      const float throughput = path_throughput(path);
      if(throughput > 0.0f && (path->v[v2].mode & s_emit))
      {
        weight = sampler_mis(path, 1);
        pointsampler_splat(path, throughput * weight / path->mmlt.p);
        _configuration_accum(rt.sampler->conf, npt, nlt, throughput * weight/path->mmlt.p);
        goto success;
      }
      goto no_path;
    }
  }

  // no light/eye path or only next event estimation? no need to connect then.
  if(npt <= 1 || nlt <= 1) goto no_path;

  // FIXME: check with kmlt, there seems to be the assumption that the light path is the first one in
  // the connected construct. probably need eye path for hslt though.
  // will append lt_path to path.
  const float throughput = path_connect(path, &lt_path);
  if(throughput <= 0.f) goto no_path; // might fail due to excessive path length
  weight = sampler_mis(path, lt_path.length);
  pointsampler_splat(path, throughput * weight / path->mmlt.p);
  _configuration_accum(rt.sampler->conf, npt, nlt, throughput * weight/path->mmlt.p);
  goto success;
no_path:
  path_init(path); // invalidate
  return;
success:
  assert(path->mmlt.p > 0.0);
  // account for additional sampling weight (so mlt will get correct mean image brightness)
  path->throughput *= weight / path->mmlt.p;
  if(path->throughput >= 0.0) return;
  else goto no_path;
}

double sampler_throughput(path_t *path)
{
  if(path->length < 2) return 0.0;
  double f = path_measurement_contribution_dx(path);
  if(f <= 0.0) return 0.0;
  // now compute sum of all possible (vertex area measure) pdfs with this technique:
  double pdf = 0.0;

  // vertex area pdf as if shot from
  float vpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  float vpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  float vpdf_nee_fwd, vpdf_nee_adj; // next event estimation pdfs on both ends

  // evaluate pdfs and adjoint pdfs (adjoint is wrt to the current path, always the opposite direction)
  vpdf_nee_fwd = path_pdf_next_event(path, path->length-1);
  vpdf_nee_adj = path_pdf_next_event_adjoint(path, 0);
  for(int k=0;k<path->length;k++)
  {
    vpdf_fwd[k] = path_pdf_extend(path, k);
    vpdf_adj[k] = path_pdf_extend_adjoint(path, k);
  }

  // need to normalise probability to sample this path configuration
  // to this particular path length:
  float conf_p_sum = 1.0f;
#if 0 // would need that, but normalises away in acceptance probability anyways.
  for(int num_fwd=0;num_fwd<=path->length;num_fwd++)
  {
    if(path->v[0].mode & s_sensor)
      conf_p_sum += _configuration_p(rt.sampler->conf, num_fwd, path->length-num_fwd);
    else
      conf_p_sum += _configuration_p(rt.sampler->conf, path->length-num_fwd, num_fwd);
  }
#endif
  assert(conf_p_sum > 0.0);

  // num_fwd is how many vertices were created from the same side as we were constructed
  for(int num_fwd=0;num_fwd<=path->length;num_fwd++)
  {
    // specular connections have zero probability (but pdf() will return non-zero)
    if((num_fwd > 0 && num_fwd < path->length) &&
       ((path->v[num_fwd-1].mode & s_specular) ||
        (path->v[num_fwd].mode & s_specular)))
      continue;

    double opdf = 1.0f;
    if(num_fwd == 0)
    {
      // all vertices created from the other side.
      if(path->v[0].mode & s_sensor) opdf = 0.0f; // lt can't hit the lens
      else for(int k=0;k<path->length-1;k++) opdf *= vpdf_adj[k]; // pt, last two created in one path_extend() call
    }
    else if(num_fwd == 1)
    {
      // next event estmiation from the other side
      for(int k=1;k<path->length-1;k++) opdf *= vpdf_adj[k]; // len-2 and len-1 created in one call again
      opdf *= vpdf_nee_adj;
    }
    else if(num_fwd == path->length-1)
    {
      // next event estimation from our side
      for(int k=1;k<num_fwd;k++) opdf *= vpdf_fwd[k]; // start at 1, first two coupled.
      opdf *= vpdf_nee_fwd;
    }
    else if(num_fwd == path->length)
    {
      // scatter and hit by chance. lt doesn't do that
      if(path->v[0].mode & s_emit) opdf = 0.0f;
      else for(int k=1;k<path->length;k++) opdf *= vpdf_fwd[k]; // pt, first two created in one call
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      for(int k=1;k<num_fwd;k++) opdf *= vpdf_fwd[k];
      for(int k=num_fwd;k<path->length-1;k++) opdf *= vpdf_adj[k];
    }
    // opdf is now pdf of other technique, we need to also take into account the proability to select this
    // strategy in this thread:
    if(path->v[0].mode & s_sensor)
      pdf += fabs(opdf)*_configuration_p(rt.sampler->conf, num_fwd, path->length-num_fwd)/conf_p_sum;
    else
      pdf += fabs(opdf)*_configuration_p(rt.sampler->conf, path->length-num_fwd, num_fwd)/conf_p_sum;
  }

  if(pdf <= 0.0) return 0.0;
  return f/pdf;
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : bi-directional path tracer with just one connection\n");
}
