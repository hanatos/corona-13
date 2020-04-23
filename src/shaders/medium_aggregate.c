
#include "pathspace.h"
#include "points.h"
#include "pointsampler.h"
#include "shader.h"
#include <assert.h>

// aggregates multiple heterogeneous volumes in one combined shader.
// implements additive collision coefficients:
// - sampling samples all and returns shortest free flight distance
// - transmittances multiply
// - pdfs are computed via transmittance * sum(mu_t)
// - phase functions are interpreted as convex combinations via combined scattering albedos ( mu_s/sum(mu_t) )

typedef struct
{
  mf_t mu_s;
  mf_t mu_t;
}
cache_entry_t;

typedef struct
{
  cache_entry_t v[PATHSPACE_MAX_VERTS];
}
path_cache_t;

typedef struct
{
  int curr;
  path_cache_t p[5];
}
cache_t;

// __thread cache_t medium_aggregate_tls;

typedef struct
{
  // TODO: use some hierarchy
  int cnt;
  int child[5];
  int mshader;
}
aggr_t;

float prepare(path_t *p, int v, void *data)
{
  aggr_t *a = data;
  // TODO: need extra pointer backdoor on path_t vertex
  // TODO: make it point to some kind of thread local memory block in here.
  // TODO: could pull out via path_t pointer.
  // TODO: support 5 paths per thread and do round robin?

  p->v[v].shading.em = mf_set1(0.0f);
  if(primid_invalid(p->v[v].hit.prim))
  {
    // volume lookup in all sub-volumes, count number
    int num = 0;
    for(int k=0;k<a->cnt;k++)
    {
      shader_so_t *child = rt.shader->shader + a->child[k];
      child->prepare(p, v, child->data);
      if(mf_any(mf_gt(p->v[v].interior.mu_t, mf_set1(0.0f))))
      {
        p->v[v].interior.l_mu_t[num] = p->v[v].interior.mu_t;
        p->v[v].interior.l_mu_s[num] = p->v[v].interior.mu_s;
        p->v[v].interior.l_g[num] = p->v[v].interior.mean_cos;
        num++;
      }
    }
    assert(num <= 5); // path space doesn't support more right now.
    p->v[v].interior.num_lobes = num;
    // now aggregate into interior
    p->v[v].interior.mu_t = mf_set1(0.0f);
    p->v[v].interior.mu_s = mf_set1(0.0f);
    p->v[v].interior.mean_cos = 0.0f;
    for(int k=0;k<num;k++) p->v[v].interior.mu_t = mf_add(p->v[v].interior.mu_t, p->v[v].interior.l_mu_t[k]);
    for(int k=0;k<num;k++) p->v[v].interior.mu_s = mf_add(p->v[v].interior.mu_s, p->v[v].interior.l_mu_s[k]);
    for(int k=0;k<num;k++)
      p->v[v].interior.mean_cos += mf(p->v[v].interior.l_mu_s[k], 0) / mf(p->v[v].interior.mu_s, 0) * p->v[v].interior.l_g[k];
  }
  p->v[v].interior.shader = a->mshader;
  return 1.0f;
}

int volume_enabled(void *data)
{
  return 1;
}

int shape_init(uint32_t shapeid, shader_so_t *self)
{
#if 1 // called on the individual volumes already, nothing to do for us
  aggr_t *a = self->data;
  a->mshader = self - rt.shader->shader;
  if(a->mshader != rt.shader->exterior_medium_shader) return 0; // we want to rescale and discard only if we're the exterior medium
  int discard = 0;
  for(int k=0;k<a->cnt;k++)
  {
    rt.shader->exterior_medium_shader = a->child[k];
    if(rt.shader->shader[a->child[k]].shape_init)
      discard |= rt.shader->shader[a->child[k]].shape_init(shapeid, rt.shader->shader + a->child[k]);
  }
  rt.shader->exterior_medium_shader = a->mshader;
  return discard;
#else
  return 0;
#endif
}

int init(FILE* f, void** data)
{
  // get indices to children
  aggr_t *a = calloc(1, sizeof(*a));
  *data = a;

  int dreggn = 0;
  // read sample shader number
  if(fscanf(f, "%d", &(a->cnt)) < 1)
  {
    fprintf(stderr, "[medium_aggregate] could not read all parameters! expecting: <num> <s0> <s1>..\n");
    return 1;
  }
  assert(a->cnt <= 5);
  assert(a->cnt > 0);
  for(int k=0;k<a->cnt;k++)
  {
    if(fscanf(f, "%d", a->child + k) < 1)
    {
      fprintf(stderr, "[medium_aggregate] failed to read %d-th shader num!\n", k);
      return 1;
    }
  }
  dreggn = fscanf(f, "%*[^\n]\n");
  shader_so_t *self = (shader_so_t *)data;
  assert(self->data == *data);

  // get our shader number and offset the others if negative numbers:
  int us = self - rt.shader->shader;
  for(int k=0;k<a->cnt;k++) if(a->child[k] < 0) a->child[k] = us + a->child[k];

  return dreggn == -1;
}

void cleanup(void *data)
{
  // free our struct, every child gets their own cleanup call already
  aggr_t *a = data;
  free(a);
}

// distance sampling and transmittance
float volume_sample(path_t *p, int e, void *data)
{
  // iterate over all volumes, return shortest distance
  const aggr_t *a = data;
  float old_dist = p->e[e].dist;
  float dist = old_dist;
  mf_t transmittance[5] = {mf_set1(0.0f)};
  int valid[5] = {0};
  mf_t mu_t[5] = {mf_set1(0.0f)};
  for(int k=0;k<a->cnt;k++)
  {
    // if using halton points, we need to fake randoms for all but the first medium (to decorrelate our estimators):
    if(k) pointsampler_enable_fake_random(rt.pointsampler);
    if(k) pointsampler_set_fake_random(rt.pointsampler, s_dim_free_path, points_rand(rt.points, common_get_threadid()));
    float new_dist = rt.shader->shader[a->child[k]].volume_sample(p, e, rt.shader->shader[a->child[k]].data);
    if(k) pointsampler_disable_fake_random(rt.pointsampler);
    if(new_dist < dist)
    {
      // we got a new distance, all other transmittances and mu_t are invalid now:
      for(int i=0;i<k;i++) valid[i] = 0;
      dist = new_dist;
      p->e[e].dist = new_dist; // clip closer, also need to set for transmittance
    }
    valid[k] = 1;
    mu_t[k] = p->v[e].interior.mu_t;
    transmittance[k] = p->e[e].transmittance;
  }
  if(dist < FLT_MAX) for(int k=0;k<a->cnt;k++)
  {
    if(!valid[k])
    { // ouch.
      shader_so_t *child = rt.shader->shader + a->child[k];
      transmittance[k] = child->volume_transmittance(p, e, child->data);
      for(int i=0;i<3;i++)
        p->v[e].hit.x[i] = p->v[e-1].hit.x[i] + p->e[e].omega[i] * dist;
      child->prepare(p, e, child->data);
      mu_t[k] = p->v[e].interior.mu_t;
    }
  }
  p->v[e].interior.mu_t = mf_set1(0.0f);
  p->e[e].transmittance = mf_set1(1.0f);
  if(dist < FLT_MAX) for(int k=0;k<a->cnt;k++)
  {
    if(mf_any(mf_lt(transmittance[k], mf_set1(1.0f))))
    {
      p->v[e].interior.mu_t = mf_add(p->v[e].interior.mu_t, mu_t[k]);
      p->e[e].transmittance = mf_mul(p->e[e].transmittance, transmittance[k]);
    }
  }
  // update_throughput in pathspace depends on transmittance / edge pdf to evaluate the contribution:
  if(dist == FLT_MAX)
    p->e[e].transmittance = mf_set1(1.0f);
  if(dist < old_dist)
    p->e[e].pdf = mf_mul(p->e[e].transmittance, p->v[e].interior.mu_t);
  else
    p->e[e].pdf = p->e[e].transmittance;
  p->e[e].dist = old_dist; // restore, we're only returning a proposition:
  return dist;
}

mf_t volume_transmittance(path_t *p, int e, void *data)
{
  const aggr_t *a = data;
  mf_t transmittance = mf_set1(1.0f);
  for(int k=0;k<a->cnt;k++)
    transmittance = mf_mul(transmittance, rt.shader->shader[a->child[k]].volume_transmittance(p, e, rt.shader->shader[a->child[k]].data));
  p->e[e].transmittance = transmittance;
  p->e[e].pdf = mf_set1(1.0f); // was deterministic
  return transmittance;
}

mf_t volume_pdf_adj(path_t *p, int e, void *data)
{
  // return transmittance * mu_t at start point
  if(mf_any(mf_gt(p->e[e].transmittance, mf_set1(0.0f)))) // assume good cached value
  {
    if(p->v[e-1].material_modes & s_volume) return mf_mul(p->e[e].transmittance, p->v[e-1].interior.mu_t);
    return p->e[e].transmittance;
  }
  float dist = p->e[e].dist;
  if(p->v[e].flags & s_environment)
  { // compute real distance between the points:
    float delta[3] = {
      p->v[e].hit.x[0] - p->v[e-1].hit.x[0],
      p->v[e].hit.x[1] - p->v[e-1].hit.x[1],
      p->v[e].hit.x[2] - p->v[e-1].hit.x[2]};
    dist = sqrtf(dotproduct(delta, delta));
  }
  const mf_t transmittance = volume_transmittance(p, e, data);
  if(p->v[e-1].material_modes & s_volume) return mf_mul(transmittance, p->v[e-1].interior.mu_t);
  return transmittance;
}

mf_t volume_pdf(path_t *p, int e, void *data)
{
  // return transmittance * mu_t at end point
  if(mf_any(mf_gt(p->e[e].transmittance, mf_set1(0.0f)))) // assume good cached value
  {
    if(p->v[e].material_modes & s_volume) return mf_mul(p->e[e].transmittance, p->v[e].interior.mu_t);
    return p->e[e].transmittance;
  }
  const mf_t transmittance = volume_transmittance(p, e, data);
  if(p->v[e].material_modes & s_volume) return mf_mul(transmittance, p->v[e].interior.mu_t);
  return transmittance;
}

// phase function
mf_t pdf(path_t *p, int e1, int v, int e2, void *data)
{
  // eval all pdf
  aggr_t *a = data;
  int num = p->v[v].interior.num_lobes;
  if(num <= 0) return mf_set1(0.0f); // could happen if mu_t == 0 is detected too late
  mf_t sum_pdf = mf_set1(0.0f);
  for(int k=0;k<num;k++)
  {
    shader_so_t *shader = rt.shader->shader + a->child[k];
    sum_pdf = mf_fma(shader->pdf(p, e1, v, e2, shader->data), mf_div(p->v[v].interior.l_mu_s[k], p->v[v].interior.mu_s), sum_pdf);
  }
  return sum_pdf;
}

mf_t brdf(path_t *p, int v, void *data)
{
  // eval all phase functions and weight by respective mu_s
  aggr_t *a = data;
  int num = p->v[v].interior.num_lobes;
  if(num <= 0) return mf_set1(0.0f); // could happen if mu_t == 0 is detected too late
  const mf_t mu_s = p->v[v].interior.mu_s;
  mf_t phase = mf_set1(0.0f);
  for(int k=0;k<num;k++)
  {
    shader_so_t *shader = rt.shader->shader + a->child[k];
    p->v[v].interior.mu_s = p->v[v].interior.l_mu_s[k]; // pretend single lobe
    phase = mf_add(phase, shader->brdf(p, v, shader->data));
  }
  p->v[v].interior.mu_s = mu_s; // restore
  return phase;
}

mf_t sample(path_t *p, void *data)
{
  aggr_t *a = data;
  // sample convex combination: mu_s / sum(mu_s)
  const float rand = pointsampler(p, s_dim_scatter_mode);
  int v = p->length - 1;
  int num = p->v[v].interior.num_lobes;
  if(num <= 0) return mf_set1(0.0f); // could happen if mu_t == 0 is detected too late
  int lobe = num-1;
  float cdf = 0.0;
  for(int k=0;k<num-1;k++)
  {
    cdf += mf(p->v[v].interior.l_mu_s[k], 0) / mf(p->v[v].interior.mu_s, 0);
    if(rand < cdf)
    {
      lobe = k;
      break;
    }
  }
  shader_so_t *shader = rt.shader->shader + a->child[lobe];
  shader->sample(p, shader->data);
  // eval pdf: save one computation by pulling out p->v[v+1].pdf after sample()
  mf_t pdf = mf_mul(p->v[v+1].pdf, p->v[v].interior.l_mu_s[lobe]);
  for(int k=0;k<num;k++)
  {
    if(k == lobe) continue;
    shader_so_t *shader = rt.shader->shader + a->child[k];
    pdf = mf_fma(shader->pdf(p, v, v, v+1, shader->data), p->v[v].interior.l_mu_s[k], pdf);
  }

  p->v[v+1].pdf = mf_div(pdf, p->v[v].interior.mu_s);
  // return w; // weights cancel with cdf above, so it's just phi_i * sum(mu_s)
  return p->v[v].interior.mu_s; // == brdf(p, v, data) / p->v[v+1].pdf;
}

// next event estimation
mf_t volume_pdf_nee(const path_t *p, int v, void *data)
{
  // need to ask all of them and sum up
	assert(0 && "not implemented");
	return mf_set1(0.f);
}

mf_t volume_sample_nee(path_t *p, int v, void *data)
{
  // choose one and sample that
	assert(0 && "not implemented");
	return mf_set1(0.f);
}

// light tracing
mf_t volume_sample_light(path_t *p, void *data)
{
	assert(0 && "not implemented");
	return mf_set1(0.f);
}

mf_t volume_pdf_light(const path_t *p, int v, void *data)
{
	assert(0 && "not implemented");
	return mf_set1(0.f);
}
