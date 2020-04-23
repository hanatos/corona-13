#pragma once

#include "halfvec.h"

static inline float multichain_perturb_distance(
    path_t *curr,
    path_t *tent,
    int v)                 // e[v+1].dist will be mutated (the outgoing distance of vertex v)
{
  float T = 1.0f;
  // assume no volume interaction:
  tent->e[v+1].dist = FLT_MAX;
  // mutate distances: mutate all but the last segment. if an edge reaches
  // the surface again (be it for total internal reflection or
  // transmission), this distance cannot be sampled here but will be
  // dictated by ray tracing.

  // FIXME: need to estimate this carefully, especially in case v is not in volume and v+1 is. in this case
  // a distance with mu from v+1 would be much shorter if chosen uniformly according to this.
  const float curr_mu_t = curr->v[v+1].interior.mu_t;
  const float sigma_step = 4.0;

  if(curr->v[v+1].mode & s_volume)
  { // mutate distance in scattering medium, asymmetric pdf!
#if 0
    // obtain hypothetical random number:
    const float curr_rand = 1.0 - expf(-curr->e[v+1].dist * curr_mu_t);
    const float tent_rand = sample_mutate_rand(curr_rand, points_rand(rt.points, common_get_threadid()), 0.01f);
    tent->e[v+1].dist = - logf(1.0f-tent_rand)/curr_mu_t;
    // may consider going through the shader_vol_sample() interface.
#endif
    float g0, g1;
    const int tid = common_get_threadid();
    const float r0 = points_rand(rt.points, tid);
    const float r1 = points_rand(rt.points, tid);
    sample_gaussian(r0, r1, &g0, &g1);
    tent->e[v+1].dist = curr->e[v+1].dist + g0 / (sigma_step*curr_mu_t);
    if(tent->e[v+1].dist < 0.0) return 0.0;
  }

  tent->v[v+1].mode = curr->v[v+1].mode;
  tent->e[v+1].transmittance = 0.0f;
  if(path_propagate(tent, v+1, s_propagate_mutate) ||
      (tent->v[v+1].flags                    != curr->v[v+1].flags) ||
      (tent->v[v+1].material_modes           != curr->v[v+1].material_modes) ||
      (tent->v[v+1].hit.shader               != curr->v[v+1].hit.shader) ||
      (tent->v[v+1].interior.shader          != curr->v[v+1].interior.shader) ||
      (primid_invalid(tent->v[v+1].hit.prim) != primid_invalid(curr->v[v+1].hit.prim))) 
    return 0.0f;

  // transition probability tent f/T  / curr f/T => curr_p / tent_p
  if(tent->v[v+1].mode & s_volume)
  {
    // fprintf(stderr, "mut vol mu_t %g %g dist %g->%g p = %g\n",
    //     curr_mu_t, tent_mu_t,
    //     curr->e[v+1].dist, tent->e[v+1].dist,
    // curr_mu_t / tent_mu_t *
    //   expf(-(
    //     curr->e[v+1].dist * curr_mu_t -
    //     tent->e[v+1].dist * tent_mu_t)));
    const float tent_mu_t = tent->v[v+1].interior.mu_t;
    float delta = curr->e[v+1].dist - tent->e[v+1].dist;
    return exp(-delta*delta/ (2.0*sigma_step*tent_mu_t*sigma_step*tent_mu_t));
#if 0
    return curr_mu_t / tent_mu_t *
      expf(-(
        curr->e[v+1].dist * curr_mu_t -
        tent->e[v+1].dist * tent_mu_t));
#endif
  }
  return T;
}

// perturb outgoing direction in solid angle domain.
// i suppose this would be better in random number domain.
static inline float multichain_perturb_dwp(
    path_t *curr,
    path_t *tent,
    int v)
{
  if(curr->v[v].mode & s_specular)
  { // fix specular normal, avoid drift:
    tent->v[v].diffgeo.h[0] = 0.0f;
    tent->v[v].diffgeo.h[1] = 0.0f;
    tent->v[v].diffgeo.h[2] = 1.0f;
    if(halfvec_reflect(tent, v)) return 0.0f;
    return 1.0f; // symmetric, stay at constraint
  }
  // special case for envmaps
  if(curr->v[v].flags & s_environment) return 1.0f; // could perturb position instead, but wtf.

  const int tid = common_get_threadid();

  if(v == 0) goto generic_solid_angle; // should only happen for light tracing

  float r_omega_x, r_omega_y, r_scatter_mode;
  if(shader_inverse_sample(curr, v, &r_omega_x, &r_omega_y, &r_scatter_mode))
    return 0.0f;
  if(r_omega_x < 0.0) // random number inversion not supported
    goto generic_solid_angle;

  // perturb random numbers:
  if(curr->v[v].mode & s_volume)
  {
    r_omega_x = sample_mutate_rand(r_omega_x, points_rand(rt.points, tid), 0.001f);
    r_omega_y = sample_mutate_rand(r_omega_y, points_rand(rt.points, tid), 0.001f);
  }
  else
  {
    r_omega_x = sample_mutate_rand(r_omega_x, points_rand(rt.points, tid), 0.01f);
    r_omega_y = sample_mutate_rand(r_omega_y, points_rand(rt.points, tid), 0.01f);
  }
  // actually i think we want to leave this as it is:
  // r_scatter_mode = sample_mutate_rand(r_scatter_mode, points_rand(rt.points, tid), 0.05f);

  // put pointsampler into fake random state for our thread:
  pointsampler_enable_fake_random(rt.pointsampler);
  pointsampler_set_fake_random(rt.pointsampler, s_dim_omega_x, r_omega_x);
  pointsampler_set_fake_random(rt.pointsampler, s_dim_omega_y, r_omega_y);
  pointsampler_set_fake_random(rt.pointsampler, s_dim_scatter_mode, r_scatter_mode);
  // sample new outgoing direction
  const int oldlen = tent->length;
  tent->length = v;
  // const float weight_tent =
  shader_sample(tent);
  tent->length = oldlen; // reset
  // reset pointsampler state:
  pointsampler_disable_fake_random(rt.pointsampler);
  const float pdf_curr = shader_pdf(curr, v);
  return pdf_curr / tent->v[v+1].pdf;

  // fall back to generic solid angle distribution
generic_solid_angle:;
  if(v == 0 && curr->v[v].flags & s_environment)
  { // perturb vertex area position, too:
    const float *aabb = accel_aabb(rt.accel);
    float a[3], b[3];
    get_onb(curr->e[1].omega, a, b);
#if 0 // do this? need to make symmetric?
    float mina = INFINITY, maxa = - INFINITY, minb = INFINITY, maxb = - INFINITY;
    for(int i=0;i<4;i+=3) for(int j=0;j<4;j+=3) for(int k=0;k<4;k+=3) 
    {
      const float dota = aabb[i+0]*a[0] + aabb[j+1]*a[1] + aabb[k+2]*a[2];
      const float dotb = aabb[i+0]*b[0] + aabb[j+1]*b[1] + aabb[k+2]*b[2];
      if(dota < mina) mina = dota;
      if(dota > maxa) maxa = dota;
      if(dotb < minb) minb = dotb;
      if(dotb > maxb) maxb = dotb;
    }
#endif
    const float far = MAX(aabb[5] - aabb[2],
                      MAX(aabb[4] - aabb[1],
                          aabb[3] - aabb[0]));
    float g0, g1;
    const int tid = common_get_threadid();
    const float r0 = points_rand(rt.points, tid);
    const float r1 = points_rand(rt.points, tid);
    sample_gaussian(r0, r1, &g0, &g1);
    for(int k=0;k<3;k++)
      tent->v[0].hit.x[k] = curr->v[0].hit.x[k]
        + g0 * 0.02*far*a[k] + g1 * 0.02*far*b[k];
  }

  // assure symmetry:
  const float g = .5f*(tent->v[v].interior.mean_cos + curr->v[v].interior.mean_cos);
  const float r = .5f*(tent->v[v].shading.roughness + curr->v[v].shading.roughness);
  float k = 0.0f;
  if(curr->v[v].mode & s_volume) // volumetric scattering:
    k = 5000.f * g; // XXX TODO: needs testing and adjustment
  else if(curr->v[v].mode & s_fiber) // hair:
    k = 20.f/(r*r + 1e-4f); // XXX TODO: needs testing and adjustment
  else if(v == 0 && (curr->v[0].mode & s_emit))
    k = 2000.f/(r*r + 1e-10f);
  else // surface case:
    k = 5000.f/(r*r + 1e-4f);

  // float wc[3], wt[3]; // directions in tangent frame
  // wc[0] = dotproduct(curr->v[v].hit.a, curr->e[v+1].omega);
  // wc[1] = dotproduct(curr->v[v].hit.b, curr->e[v+1].omega);
  // wc[2] = dotproduct(curr->v[v].hit.n, curr->e[v+1].omega);
  float a[3], b[3], x, y, z;
  // get_onb(wc, a, b);
  get_onb(curr->e[v+1].omega, a, b);
  sample_cos_k(&x, &y, &z, k,
      points_rand(rt.points, tid), points_rand(rt.points, tid));
  for(int i=0;i<3;i++)
    tent->e[v+1].omega[i] = x * a[i] + y * b[i] + z * curr->e[v+1].omega[i];
  // for(int i=0;i<3;i++)
  //   wt[i] = x * a[i] + y * b[i] + z * wc[i];
  // for(int i=0;i<3;i++)
  //   tent->e[v+1].omega[i] = wt[0] * tent->v[v].hit.a[i]
  //     + wt[1] * tent->v[v].hit.b[i]
  //     + wt[2] * tent->v[v].hit.n[i];
  normalise(tent->e[v+1].omega);

  // return ratio of outgoing cosines (or sin for hair or nothing in volumes..) to be in dwp:
  return path_lambert(tent, v, tent->e[v+1].omega) / path_lambert(curr, v, curr->e[v+1].omega);
}

// returns dx acceptance probability f/T / f/T for this sub-path.
static inline double multichain_perturb(
    path_t *curr,           // current path
    path_t *tent,           // this path will be modified, needs to have inited tent->sensor struct
    const int begin,        // first vertex (e.g. 0 for start on camera)
    const int end)          // last moving vertex, end of chain (end + 1 stays fixed)
{
  if(end == begin) return 1.0;
  double T = 1.0;
  // retrace segment (a..b) by lens perturbation + multichain half vector perturbation
  // trace tent path keeping half vectors to re-initialise
  // all hit point infos, medium transactions and to get new error vector

  // re-evaluate shading parameters on start point in case wavelength has changed on the outside
  if(begin == 0)
  {
    if(curr->v[0].mode & s_sensor)
    {
      shader_exterior_medium(tent);
      // sensor modified from outside
      if(!(tent->sensor.pixel_i >= 0 && tent->sensor.pixel_j >= 0 &&
           tent->sensor.pixel_i < view_width() && tent->sensor.pixel_j < view_height()))
        return 0.0f;

      // sample new outgoing direction and corresponding dwp pdf.
      // G required to convert p(w_1) to p(x_1) cancels with measurement.
      // we're actually interested in the ratio of perturbation pdf, but
      // these are composed of p() * J and p() is symmetric.
      float throughput = view_cam_sample(tent);
      if(throughput <= 0.0) return 0.0f;
      T *= view_cam_pdf(curr, 0) / tent->v[1].pdf;
    }
    else
    {
      shader_prepare(tent, begin);
      shader_exterior_medium(tent);
      // mutate outgoing direction on light source, account for change in measurement:
      T *= multichain_perturb_dwp(curr, tent, 0);
      if(!(T > 0)) return 0.0;
    }
  }
  else
  {
    shader_prepare(tent, begin);
    T *= multichain_perturb_dwp(curr, tent, begin);
    if(!(T > 0)) return 0.0;
    T *= shader_brdf(tent, begin) / shader_brdf(curr, begin);
    if(!(T > 0)) return 0.0;
  }

  T *= multichain_perturb_distance(curr, tent, 0);
  if(!(T > 0)) return 0.0;

  // create segment (a..b) by tracing with slightly perturbed half vectors.
  // loop over vertices in current path. will write to tent->v[v+1].
  for(int v=1;v<end;v++)
  {
    T *= multichain_perturb_dwp(curr, tent, v);
    if(!(T > 0)) return 0.0;
    T *= multichain_perturb_distance(curr, tent, v);
    if(!(T > 0)) return 0.0;
  }
  return T * path_measurement_contribution_dwp(tent, begin, end) /
             path_measurement_contribution_dwp(curr, begin, end);
}

// same as perturb, but also connect to the fixed vertex end+1.
// this implements an extended lens perturbation.
static inline double multichain_perturb_connect(
    path_t *curr,           // current path
    path_t *tent,           // this path will be modified, needs to have inited tent->sensor struct
    const int end)          // last moving vertex, end of chain (end + 1 stays fixed)
{
  assert(end <= curr->length);

  // TODO: optimise this copy away! apparently we need to set shader id somewhere:
  *tent = *curr;
  const int tid = common_get_threadid();
  float g1, g2;
  const float r1 = points_rand(rt.points, tid);
  const float r2 = points_rand(rt.points, tid);
  sample_gaussian(r1, r2, &g1, &g2);
  const float px = 10.f; // this is one sigma width of the jump
  tent->sensor.pixel_i += g1 * px;
  tent->sensor.pixel_j += g2 * px;
  // mutate point on aperture, after mutating outgoing direction (pixel)
  tent->sensor.aperture_x = sample_mutate_rand(tent->sensor.aperture_x, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_y = sample_mutate_rand(tent->sensor.aperture_y, points_rand(rt.points, tid), 0.3f);
  tent->sensor.aperture_set = 1;
  tent->sensor.pixel_set = 1;
  tent->lambda = curr->lambda;
  // XXX would need to eval the rest of the path, too!
  // tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);
  tent->time = curr->time; // <= TODO

  double T = multichain_perturb(curr, tent, 0, end);
  if(end == tent->length) return T; // no connection needed
  if(!(T > 0.0)) return 0.0;
  
  // test connection
  tent->v[end+1].mode = curr->v[end+1].mode;
  tent->e[end+1].transmittance = 0.0f;
  if(path_project(tent, end+1, s_propagate_mutate) ||
      (tent->v[end+1].flags           != curr->v[end+1].flags) ||
      (tent->v[end+1].hit.shader      != curr->v[end+1].hit.shader) ||
      (tent->v[end+1].interior.shader != curr->v[end+1].interior.shader) ||
      (primid_invalid(tent->v[end+1].hit.prim) != primid_invalid(curr->v[end+1].hit.prim))) 
    return 0.0f;

  // check whether we actually arrived at vertex c
  for(int k=0;k<3;k++)
    if(fabsf(tent->v[end+1].hit.x[k] - curr->v[end+1].hit.x[k]) > HALFVEC_REL_SPATIAL_EPS *
        MAX(MAX(fabsf(tent->v[end+1].hit.x[k]), fabsf(curr->v[end+1].hit.x[k])), 1.0))
      return 0.0f;
  T *= shader_brdf(tent, end) * shader_vol_transmittance(tent, end+1) * path_G(tent, end+1);
  T /= shader_brdf(curr, end) * shader_vol_transmittance(curr, end+1) * path_G(curr, end+1);
  // T *= shader_vol_transmittance(tent, end+1);// * path_G(tent, end+1);
  // T /= shader_vol_transmittance(curr, end+1);// * path_G(curr, end+1);
  // return T;// XXX the below seems to not work (always rejects)
  shader_prepare(tent, end+1); // potentially update due to direction/wavelength change
  if(end+1 == tent->length-1)
  {
    T *= lights_eval_vertex(tent, tent->length-1);
    T /= lights_eval_vertex(curr, curr->length-1);
  }
  else
  {
    T *= shader_brdf(tent, end+1);
    T /= shader_brdf(curr, end+1);
  }
  // nan-check:
  if(!(T >= 0.0)) return 0.0;
  return T;
}
