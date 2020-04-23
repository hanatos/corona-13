#pragma once

#include "points.h"
#include "pathspace.h"
#include "pathspace/nee.h"
#include "pathspace/halfvec.h"
#include "pathspace/multichain.h"
#include "shader.h"
#include "camera.h"

typedef struct vmlt_mmlt_t
{
}
vmlt_mmlt_t;

void *mmlt_init()
{
  vmlt_mmlt_t *d = malloc(sizeof(*d));
  memset(d, 0, sizeof(*d));
  return d;
}

void mmlt_cleanup(void *data)
{
  vmlt_mmlt_t *d = data;
  free(d);
}

float mmlt_suitability(const path_t *p, void *data)
{
  // can't construct path from scratch, so demand initialized sample (length > 0)
  if(p->length < 2) return 0.0f;
  return 1.0f;
}

static inline void _mmlt_c_cdf(const path_t *p, int b, float *cdf)
{
  // don't want to place c at vert 1, ever.
  b = MAX(1, b);
  for(int k=0;k<=b;k++) cdf[k] = 0.0f;
  // get c ~ roughness
  for(int k=b+1;k<p->length-1;k++)
  {
    if(p->v[k].mode & s_fiber)
    { // don't want a half vec chain to span across fiber scattering
      // events, we don't have good derivatives for it (nor do i
      // expect visibility of such connections to play nicely)
      cdf[k++] = 1e10f;
      for(;k<p->length;k++) cdf[k] = 0.0f;
      break;
    }
    int indexmatched = fabsf(path_eta_ratio(p, k) - 1.0f) < 1e-3f;
    float step = p->v[k].shading.roughness;
    const float g = p->v[k].interior.mean_cos;
    if(p->v[k].mode & s_volume) step = sqrtf(1.0f - g*g);
    cdf[k] = indexmatched ? 0 : step;
  }
  cdf[p->length-1] = p->v[p->length-1].shading.roughness;
  if(b+1 < p->length)
    cdf[b+1] *= 1.5f;

  // make cdf:
  for(int k=b+2;k<p->length;k++)
    cdf[k] += cdf[k-1];
  for(int k=b+1;k<p->length-1;k++)
    cdf[k] /= cdf[p->length-1];
  cdf[p->length-1] = 1.0f;
}

static inline int _mmlt_breakup(path_t *tent, float *cdf)
{
  // pdfs are touchy, need to compute it in double and write back to float only
  // after normalisation:
  double cdfd[PATHSPACE_MAX_VERTS] = {0.0};
  // want to explore two-vertex paths with full multichain:
  if(tent->length < 3)
  {
    cdf[0] = 0.0f; cdf[1] = 1.0f;
    return 0;
  }
  // compute bidirectional pdfs:
  double tpdf_fwd[PATHSPACE_MAX_VERTS]; // our way
  double tpdf_adj[PATHSPACE_MAX_VERTS]; // opposite direction
  const float vpdf_nee_fwd = nee_pdf(tent, tent->length-1);
  const float vpdf_nee_adj = nee_pdf_adjoint(tent, 0);
  for(int k=0;k<tent->length;k++)
  {
    tpdf_fwd[k] = path_pdf_extend(tent, k);
    tpdf_adj[k] = path_pdf_extend_adjoint(tent, k);
  }
  assert(tpdf_fwd[0] == 1.0);
  assert(tpdf_adj[tent->length-1] == 1.0);
  for(int k=1;k<tent->length;k++)    tpdf_fwd[k] *= tpdf_fwd[k-1];
  for(int k=tent->length-2;k>=0;k--) tpdf_adj[k] *= tpdf_adj[k+1];

  // get breakup point by mis weight, higher pdf is better.
  const int light_dir = tent->v[0].mode & s_emit;
  assert(!light_dir); // i have no trust this would work
  for(int k=0;k<=tent->length;k++)
  { // iterate through k=number of eye path vertices
    if((k > 0 && k < tent->length) &&
      ((tent->v[k-1].mode & s_specular) || (tent->v[k].mode & s_specular)))
        continue; // no connection possible, no contribution to pdf
    double opdf = 1.0;
    if(k == 0)
    {
      // all vertices created from the other side.
      if(!light_dir) opdf = 0.0f; // lt can't hit the lens
      else opdf = tpdf_adj[0]; // pt pdf
    }
    else if(k == 1)
    { // next event estmiation from the other side
      opdf = tpdf_adj[1] * vpdf_nee_adj;
      if(tent->length < 3) opdf = 0.0; // we don't do this
    }
    else if(k == tent->length-1)
    {
      // next event estimation from our side
      opdf = tpdf_fwd[k-1] * vpdf_nee_fwd;
      if(tent->length < 3) opdf = 0.0; // we don't do this
    }
    else if(k == tent->length)
    { // scatter and hit by chance. lt doesn't do that
      if(light_dir) opdf = 0.0f;
      else opdf = tpdf_fwd[k-1]; // pt pdf
    }
    else // regular case, the two beginnings were a single path_extend() call, respectively
    {
      opdf = tpdf_fwd[k-1] * tpdf_adj[k];
    }
    cdfd[k-1] = opdf;
    assert(opdf == opdf);
  }
  // make cdf:
  for(int k=1;k<tent->length;k++)
    cdfd[k] += cdfd[k-1];
  if(!(cdfd[tent->length-1] > 0.0)) return 1;
  assert(cdfd[tent->length-1] == cdfd[tent->length-1]);
  assert(cdfd[tent->length-1] > 0.0);
  for(int k=0;k<tent->length;k++)
    cdfd[k] /= cdfd[tent->length-1];
  for(int k=0;k<tent->length-1;k++)
    cdf[k] = cdfd[k];
  cdf[tent->length-1] = 1.0f;
  return 0;
}

float mmlt_mutate(
    path_t *curr,       // constant current sample
    path_t *tent,       // tentative sample
    void *data)
{
  // vmlt_mmlt_t *d = data;
  const int tid = common_get_threadid();
  double T = 1.0; // track transition probability along the way

  // TODO: optimise this copy away! apparently we need to set shader id somewhere:
  *tent = *curr;

  float cdf_b[PATHSPACE_MAX_VERTS];
  float cdf_c[PATHSPACE_MAX_VERTS];
  if(_mmlt_breakup(tent, cdf_b)) return 0.0f;
  int breakup = sample_cdf(cdf_b, tent->length, points_rand(rt.points, tid));
  float pdf_b_tent = breakup ? cdf_b[breakup] - cdf_b[breakup-1] : cdf_b[0];
  // if(tent->length > 2) breakup = 0; // XXX DEBUG
  // fprintf(stderr, "breakup = %d / %d\n", breakup, tent->length);
  float pdf_c_tent = 1.0f;
  int c = tent->length - 1;
#if 1
  if(breakup < tent->length - 1)
  {
    _mmlt_c_cdf(tent, breakup, cdf_c);
    c = breakup + 1 + sample_cdf(cdf_c+breakup+1, tent->length-breakup-1, points_rand(rt.points, tid));
    pdf_c_tent = c ? cdf_c[c] - cdf_c[c-1] : cdf_c[0];
  }
  assert(!(tent->v[breakup].material_modes & s_specular));
  if(breakup + 1 < curr->length)
    assert(!(tent->v[breakup+1].material_modes & s_specular));
  if(c + 1 < curr->length)
    assert(!(tent->v[c+1].material_modes & s_specular));
#endif
  // breakup = MAX(0, MIN(tent->length - 1, 0));

  // mutate wavelength and time a bit. this is done in halfvec_perturb if b==0.
  tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
  tent->time = curr->time;
  tent->sensor.pixel_set = 1;
  tent->sensor.aperture_set = 1;

  if(breakup > 0)
  { // mutate eye subpath
    tent->length = breakup + 1;
    float g1, g2;
    const float r1 = points_rand(rt.points, tid);
    const float r2 = points_rand(rt.points, tid);
    sample_gaussian(r1, r2, &g1, &g2);
    const float stepsize = .5f*HALFVEC_MUTATION_STEP;
    tent->sensor.pixel_i += stepsize * g1;
    tent->sensor.pixel_j += stepsize * g2;
    // mutate point on aperture, after mutating outgoing direction (pixel)
    tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.01f);
    tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.01f);

    // camera_mutate_aperture(rt.cam, tent, points_rand(rt.points, tid), points_rand(rt.points, tid), 0.2f);
    // T /= camera_pdf_mutate_aperture(rt.cam, curr, tent, 0.2f);
    // T *= camera_pdf_mutate_aperture(rt.cam, tent, curr, 0.2f);
    T *= multichain_perturb(curr, tent, 0, breakup);
    if(!(T > 0.0)) goto fail;
  }

  if(breakup < c-1)
  {
    // mutate light subpath
    path_t curr_lp;
    path_reverse(&curr_lp, curr);
    curr_lp.length = curr->length - breakup - 1;
    path_t tent_lp = curr_lp;
    tent_lp.lambda = tent->lambda; // use what we set on the eye path earlier on
    tent_lp.time = tent->time;
    T *= multichain_perturb(&curr_lp, &tent_lp, curr->length-1-c, curr_lp.length-1);
    if(!(T > 0.0)) goto fail;

    if(breakup == 0)
    { // next event estimation connecting to eye
      tent_lp.length = curr->length;
      tent_lp.v[tent_lp.length-1].mode = s_sensor;
      tent_lp.e[tent_lp.length-1].transmittance = 0.0f;
      if(path_project(&tent_lp, tent_lp.length-1, s_propagate_mutate)) goto fail;
      path_reverse(tent, &tent_lp);
      T *= view_cam_eval(tent)/view_cam_eval(curr);
      if(!(T > 0.0)) goto fail;
      // fprintf(stderr, "acceptance %g\n", T);
      // path_print(curr, stdout);
      // path_print(tent, stdout);
    }
    else
    { // generic deterministic connect
      if(!(path_connect(tent, &tent_lp) > 0.0f)) goto fail;
      T *= shader_brdf(tent, breakup)/shader_brdf(curr, breakup);
      if(!(T > 0.0)) goto fail;
      assert(tent->v[tent->length-1].mode & s_emit);
    }
    // ratio of deterministic connection part of the measurement contribution
    T *= path_G(tent, breakup+1) * shader_vol_transmittance(tent, breakup+1) * shader_brdf(tent, breakup+1) /
        (path_G(curr, breakup+1) * shader_vol_transmittance(curr, breakup+1) * shader_brdf(curr, breakup+1));
    if(!(T > 0.0)) goto fail;
  }
  else if(breakup == c-1)
  { // next event estimation connecting to light source
    tent->length = curr->length; // reset to full length again
    tent->v[c].mode = curr->v[c].mode;
    tent->e[c].transmittance = 0.0f;
    if(path_project(tent, c, s_propagate_mutate)) goto fail; // trace visibility and initialise edf correctly
    T *= shader_brdf(tent, breakup) * path_G(tent, breakup+1) * shader_vol_transmittance(tent, breakup+1)/
        (shader_brdf(curr, breakup) * path_G(curr, breakup+1) * shader_vol_transmittance(curr, breakup+1));
    if(!(T > 0.0)) goto fail;
    if(c < curr->length-1)
      T *= shader_brdf(tent, c) / shader_brdf(curr, c);
    else // light vertex at the end:
      T *= lights_eval_vertex(tent, c) / lights_eval_vertex(curr, c);
    if(!(T > 0.0)) goto fail;
  }
  else
  {
    // breakup == c: full multichain from eye and we're done already
    assert(c == curr->length-1);
    T *= lights_eval_vertex(tent, c) / lights_eval_vertex(curr, c);
    if(!(T > 0.0)) goto fail;
  }

  // adjust for wavelength changes: re-eval measurement in constant segment:
  if(c < tent->length-1)
  {
    for(int v=c+1;v<tent->length;v++) shader_prepare(tent, c);
    T *= path_measurement_contribution_dwp(tent, c, tent->length-1)/
         path_measurement_contribution_dwp(curr, c, curr->length-1);
    if(!(T > 0.0)) goto fail;
  }

  // TODO: compute reverse pdf the same way (pdf fwd, pdf rev, max heuristic, pdf mutations)
  // const int breakup2 = _mmlt_breakup(tent, cdf_b);
  // if(breakup2 != breakup) return 0.0; // XXX may actually have to sample this by power heuristic instead.

  float pdf_c_curr = 1.0f;
  float pdf_b_curr = 1.0f;
#if 1
  if(_mmlt_breakup(tent, cdf_b)) return 0.0f;
  pdf_b_curr = breakup ? cdf_b[breakup] - cdf_b[breakup-1] : cdf_b[0];
  if(tent->length > 2)
  {
    _mmlt_c_cdf(tent, breakup, cdf_c);
    pdf_c_curr = c ? cdf_c[c] - cdf_c[c-1] : cdf_c[0];
  }
#endif

  // paranoid volume stack checks:
  if(!(path_measurement_contribution_dwp(tent, 0, tent->length-1) > 0))
    return 0.0f;
  for(int v=1;v<tent->length-2;v++)
    if(tent->v[v].mode != curr->v[v].mode)
      return 0.0f;

  if(T * pdf_c_curr / pdf_c_tent * pdf_b_curr / pdf_b_tent > 0.0)
  {
    for(int v=1;v<tent->length-2;v++)
      if(tent->v[v].mode & s_transmit)
        assert(tent->e[v].vol.ior != tent->e[v+1].vol.ior);
    assert(tent->v[tent->length-1].mode & s_emit);
    assert(path_measurement_contribution_dwp(tent, 0, tent->length-1) > 0);
  }
  // return acceptance as f->/T->  /  f<-/T<-
  return T * pdf_c_curr / pdf_c_tent * pdf_b_curr / pdf_b_tent;
fail:
  tent->length = 0;
  return 0.0f;
}

void mmlt_print_info(FILE *f, void *data)
{
  //vmlt_mmlt_t *d = data;
  fprintf(f, "         : multiplexed metropolis in path space\n");
}
