#pragma once
#include "vol/interpolation.h"
#include "vol/shaders.h"
#include "mf.h"
#include <float.h>

// quadrature rule version, unbiased enumeration of all pierced voxels
static inline mf_t vol_trace_transmittance(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,                // distance to trace to
    const mf_t sigma_t,                  // cross section, will be multiplied to opt. thick.
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
    float time,                          // time for motion blur access
    const mf_t lambda,                   // wavelengths
    mf_t *emission,                      // return accumulated volume emission
    float *last_density)
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
#define VOL_TRACE_INIT float thickness = 0.0f;
#define VOL_TRACE_LEAF \
  const float dist = MAX(0.0f, MIN(tmax, max_dist) - MAX(0.0f, tmin));\
  thickness += result[0] * dist;\
  if(emission && result[1] > 0.0f)\
    *emission = mf_fma(mf_set1(dist*result[0]), tree->shader(0, result[1], lambda), *emission);\
  if(last_density && dist > 0.0f) *last_density = result[0];
#define VOL_TRACE_MISS return mf_set1(1.0f);
#define VOL_TRACE_RETURN return mf_exp(mf_mul(mf_set1(-thickness), sigma_t));
#include "vol/trace_impl.inc"
}

#if 0 // cannot work with hero wavelengths:
// residual ratio tracking, returns noisy estimate.
static inline mf_t vol_trace_transmittance_residual_ratio_tracking(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,                // distance to track to
    const mf_t sigma_t,
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
    float time,                          // time for motion blur access
    const mf_t lambda,                   // wavelengths
    mf_t *emission,                      // return accumulated volume emission
    float *last_density)
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  lod += 1;
  float weight = 1.0f;
  float c_thickness = 0.0f; // optical thickness of control part, collect that via regular tracking
  // do regular tracking on coarser voxel grid and sample residual ratios after each tentative step
#define VOL_TRACE_INIT float thickness = 0.0f;\
  float rand = points_rand(rt.points, common_get_threadid());\
  float log_rand = -logf(1.0f-rand)/mf(sigma_t, 0);

  // in a leaf, do regular tracking up to collision. then re-init randoms and weight transmittance.
  // use ratio tracking to estimate the residual part mu_t - mu_t_min
#define VOL_TRACE_LEAF \
  float d_max = result[0];\
  float d_min = result[1];\
  float new_thickness = thickness + (d_max-d_min) * MAX(0.0f, tmax-MAX(0.0f, tmin));\
  while(log_rand >= thickness && log_rand < new_thickness)\
  {\
    const float c0 = 1.0f-expf(-thickness*mf(sigma_t,0)), c1 = 1.0f-expf(-new_thickness*mf(sigma_t,0));\
    const float r2 = (rand - c0)/(c1-c0);\
    const float M = (c1-c0)/(1.0-c0);\
    const float dist = - logf(1.0f-r2*M)/(mf(sigma_t,0)*d_max);\
    if(tmin + dist > max_dist) return mf_mul(mf_exp(mf_mul(mf_set1(-(c_thickness + dist * d_min)), sigma_t)), mf_set1(weight));\
    const vol_payload_t *payload_data = vol_node_get_payload(tree, node[sp], idx.idx);\
    vol_index_t cidx = {0};\
    cidx.i = (int)((pos2[0] + (tmin+dist)*dir[0] - aabb[sp][0])/vwd[sp+1])&7;\
    cidx.j = (int)((pos2[1] + (tmin+dist)*dir[1] - aabb[sp][1])/vwd[sp+1])&7;\
    cidx.k = (int)((pos2[2] + (tmin+dist)*dir[2] - aabb[sp][2])/vwd[sp+1])&7;\
    const int cs = vol_node_child_static(node[sp], idx.idx);\
    vol_payload_get(payload_data, cidx.idx, time, cs, result);\
    const float d = result[0];\
    assert(d <= d_max);\
    assert(d >= d_min);\
    if(d_max > d_min)\
    {\
    const float rr = (d-d_min) / (d_max-d_min);\
    if(weight < 0.0f)\
    {\
      if(points_rand(rt.points, common_get_threadid()) < rr) return mf_set1(0.0f);\
    }\
    else\
      weight *= 1.0f - rr;\
    }\
    if(weight <= 0.0f) return mf_set1(0.0f);\
    rand = points_rand(rt.points, common_get_threadid());\
    log_rand = -logf(1.0f-rand)/mf(sigma_t, 0);\
    thickness = 0;\
    new_thickness = (d_max-d_min) * MAX(0.0, tmax-MAX(0.0, tmin+dist));\
  }\
  thickness = new_thickness;\
  c_thickness += d_min * MAX(0.0f, tmax-MAX(0.0, tmin));

#define VOL_TRACE_MISS return mf_set1(1.0f);
#define VOL_TRACE_RETURN return mf_mul(mf_exp(mf_mul(mf_set1(-c_thickness), sigma_t)), mf_set1(weight));
#include "vol/trace_impl.inc"
}

// ratio tracking, returns noisy estimate only:
static inline mf_t vol_trace_transmittance_ratio_tracking(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,                // distance to track to
    const mf_t sigma_t,
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
    float time,                          // time for motion blur access
    const mf_t lambda,                   // wavelength
    mf_t *emission,                      // return accumulated volume emission
    float *last_density)
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  lod += 1;
  mf_t transmittance = mf_set1(1.0f);
  // do regular tracking on coarser voxel grid and sample residual ratios after each tentative step
#define VOL_TRACE_INIT float thickness = 0.0f;\
  const float sig_t = mf(sigma_t, 0);\
  float rand = points_rand(rt.points, common_get_threadid());\
  float log_rand = -logf(1.0f-rand)/sig_t;

  // in a leaf, do regular tracking up to collision. then re-init randoms and weight transmittance.
#define VOL_TRACE_LEAF \
  float d_max = result[0];\
  float new_thickness = thickness + d_max * MAX(0.0f, tmax-MAX(0.0f, tmin));\
  while(log_rand >= thickness && log_rand < new_thickness)\
  {\
    const float c0 = 1.0f-expf(-thickness*sig_t), c1 = 1.0f-expf(-new_thickness*sig_t);\
    const float r2 = (rand - c0)/(c1-c0);\
    const float M = (c1-c0)/(1.0-c0);\
    const float dist = - logf(1.0f-r2*M)/(sig_t*d_max);\
    if(tmin + dist > max_dist) return transmittance;\
    /* track length estimator: */\
    /*if(tmin + dist > max_dist) return 1.0f;*/\
    const vol_payload_t *payload_data = vol_node_get_payload(tree, node[sp], idx.idx);\
    vol_index_t cidx = {0};\
    cidx.i = (int)((pos2[0] + (tmin+dist)*dir[0] - aabb[sp][0])/vwd[sp+1])&7;\
    cidx.j = (int)((pos2[1] + (tmin+dist)*dir[1] - aabb[sp][1])/vwd[sp+1])&7;\
    cidx.k = (int)((pos2[2] + (tmin+dist)*dir[2] - aabb[sp][2])/vwd[sp+1])&7;\
    const int cs = vol_node_child_static(node[sp], idx.idx);\
    vol_payload_get(payload_data, cidx.idx, time, cs, result);\
    const float d = result[0];\
    assert(d <= d_max);\
    /* track length estimator: */\
    /*if(points_rand(rt.points, common_get_threadid()) < d / d_max) return 0.0f;*/\
    transmittance = mf_mul(transmittance, mf_set1(1.0f - d / d_max));\
    if(mf_all(mf_lt(transmittance, mf_set1(0.0f)))) return mf_set1(0.0f);\
    rand = points_rand(rt.points, common_get_threadid());\
    log_rand = -logf(1.0f-rand)/sig_t;\
    thickness = 0;\
    new_thickness = d_max * MAX(0.0, tmax-MAX(0.0, tmin+dist));\
  }\
  thickness = new_thickness;

#define VOL_TRACE_MISS return mf_set1(1.0f);
#define VOL_TRACE_RETURN return transmittance;
#include "vol/trace_impl.inc"
}

// naive free-flight track-length estimator, returns very noisy estimate indeed:
static inline mf_t vol_trace_transmittance_track_length(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,                // distance to track to
    const mf_t sigma_t,
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
    float time,                          // time for motion blur access
    const mf_t lambda,                   // wavelength
    mf_t *emission,                      // return accumulated volume emission
    float *last_density)
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  lod += 1;
#define VOL_TRACE_INIT float thickness = 0.0f;\
  const mf_t mask = mf_set1(0.0f);\
  float rand = points_rand(rt.points, common_get_threadid());\
  mf_t log_rand = mf_div(mf_set1(-logf(1.0f-rand)), sigma_t);

  // in a leaf, do regular tracking up to collision. then re-init randoms and weight transmittance.
#define VOL_TRACE_LEAF \
  float d_max = result[0];\
  float new_thickness = thickness + d_max * MAX(0.0f, tmax-MAX(0.0f, tmin));\
  while(mf_any(mf_and(mf_gte(log_rand, mf_set1(thickness)),\
                      mf_lt (log_rand, mf_set1(new_thickness)))))\
  {\
    const float c0 = 1.0f-expf(-thickness*mf(sigma_t,0)), c1 = 1.0f-expf(-new_thickness*mf(sigma_t,0));\
    const float r2 = (rand - c0)/(c1-c0);\
    const float M = (c1-c0)/(1.0-c0);\
    const float dist = - logf(1.0f-r2*M)/(sigma_t*d_max);\
    if(tmin + dist > max_dist) return mf_set1(1.0f);\
    const vol_payload_t *payload_data = vol_node_get_payload(tree, node[sp], idx.idx);\
    vol_index_t cidx = {0};\
    cidx.i = (int)((pos2[0] + (tmin+dist)*dir[0] - aabb[sp][0])/vwd[sp+1])&7;\
    cidx.j = (int)((pos2[1] + (tmin+dist)*dir[1] - aabb[sp][1])/vwd[sp+1])&7;\
    cidx.k = (int)((pos2[2] + (tmin+dist)*dir[2] - aabb[sp][2])/vwd[sp+1])&7;\
    const int cs = vol_node_child_static(node[sp], idx.idx);\
    vol_payload_get(payload_data, cidx.idx, time, cs, result);\
    const float d = result[0];\
    if(points_rand(rt.points, common_get_threadid()) < d / d_max) return mf_set1(0.0f);\
    thickness = 0;\
    new_thickness = d_max * MAX(0.0, tmax-MAX(0.0, tmin+dist));\
  }\
  thickness = new_thickness;

#define VOL_TRACE_MISS return mf_set1(1.0f);
#define VOL_TRACE_RETURN return mf_select(mf_set1(0.0f), mf_set1(1.0f), mask); // XXX
#include "vol/trace_impl.inc"
}
#endif

// quadrature rule version, unbiased enumeration of all pierced voxels
static inline mf_t vol_trace_transmittance_fnee(
		const vol_tree_t *tree,
		const float *pos_ws,
		const float *dir_ws,
		const float max_dist,                // maximum distance to trace to
		const mf_t mu_t,
		vol_interpolation_t interpolation,   // access with interpolation or motion blur?
		int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
		float time,                          // time for motion blur access
		const mf_t lambda,                   // wavelength
    mf_t *emission,                      // return accumulated volume emission
    float *last_density)
{
	float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
	float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
	vol_transform_w2o(tree, pos, 1);
	vol_transform_w2o(tree, dir, 0);
#define VOL_TRACE_INIT mf_t transmittance = mf_set1(1.f); mf_t last_emission = mf_set1(0.f); if(emission) *emission = mf_set1(0.f);
#define VOL_TRACE_LEAF \
	const float dist = MAX(0.0f, MIN(tmax, max_dist) - tmin);\
	mf_t transmittance_voxel = mf_exp(mf_mul(mf_set1(-result[0]*dist), mu_t));\
  last_emission = (emission && result[0] > 0.f && result[1] > 0.f) ? mf_mul(mf_set1(result[0]), tree->shader(0, result[1], lambda)) : mf_set1(0.f);\
  if(last_density && dist > 0.0f) *last_density = result[0];\
	if(emission && result[0] > 0)\
		*emission = mf_add(*emission, mf_div(mf_mul(transmittance, last_emission), mf_mul(mf_mul(mu_t, mf_set1(result[0])), mf_sub(mf_set1(1.0f), transmittance_voxel))));\
	transmittance = mf_mul(transmittance, transmittance_voxel);
#define VOL_TRACE_MISS return mf_set1(1.0f);
#define VOL_TRACE_RETURN \
	return transmittance;
#include "vol/trace_impl.inc"
}

static inline float _vol_trace_leaf(
    const vol_tree_t *tree,
    const mf_t mu_t,
    const float log_rand,
    const float rand,
    float *thickness,
    const float density,
    const float temperature,
    const float tmin,
    const float tmax,
    const mf_t lambda,
    mf_t  *emission,
    mf_t  *transmittance,
    float *last_density,
    mf_t  *last_emission)
{
  float new_thickness = *thickness + density * MAX(0.0f, tmax-MAX(0.0f, tmin));
  if(log_rand >= *thickness && log_rand < new_thickness)
  {
    float dist = 0;
    if(new_thickness > *thickness)
    {
      const float c0 = 1.0f-expf(-*thickness*mf(mu_t,0)), c1 = 1.0f-expf(-new_thickness*mf(mu_t,0));
      const float r2 = (rand - c0)/(c1-c0);
      const float M = (c1-c0)/(1.0-c0);
      dist = - logf(1.0f-r2*M)/(mf(mu_t,0)*density);
      if(emission && temperature > 0.0f)
      {
        const mf_t Le = mf_mul(mf_set1(density), tree->shader(0, temperature, lambda));
        *emission = mf_add(*emission, mf_mul(mf_set1(dist), Le));
      }
    }
    if(last_density)  *last_density = density;
    if(last_emission) *last_emission = mf_mul(mf_set1(density), tree->shader(0, temperature, lambda));
    if(transmittance) *transmittance = mf_exp(mf_mul(mf_set1(-(*thickness + density * dist)), mu_t));
    if(emission) *emission = mf_mul(*emission, mf_mul(*transmittance, mf_mul(mf_set1(density), mu_t)));
    return MAX(1e-15f, tmin + dist);
  }
  if(emission && temperature > 0.0f)
  {
    const mf_t Le = mf_mul(mf_set1(density), tree->shader(0, temperature, lambda));
    *emission = mf_add(*emission, mf_mul(mf_set1(MAX(0.0f, tmax-MAX(0.0f, tmin))), Le));
  }
  *thickness = new_thickness;
  return -1;
}

// pdf is transmittance (*mu_t at end vertex when ending in the volume)
static inline float vol_trace_sample(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,                // maximum distance to trace to
    const mf_t mu_t,
    const float rand,                    // random number to sample distance
    vol_interpolation_t interpolation,   // access with interpolation or motion blur?
    int lod,                             // LOD: 0 is finest, 1 is 8x coarser, ..
    float time,                          // time for motion blur access
    const mf_t lambda,                   // wavelength for black body
    mf_t *emission,                     // return accumulated emission, 0 if not requested
    mf_t *transmittance,                // return transmittance if not 0
    mf_t *last_emission,                // return temperature in last voxel, if not 0
    float *last_density)                 // return density in last voxel, if not 0
{
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  // very slight speedup gained by optimising out a few exp by doing
  // cdf in log space (optical thickness).
  // M = 1-tau_voxel(tmax-tmin)
#define VOL_TRACE_INIT float thickness = 0.0f;\
  const float log_rand = -logf(1.0f-rand)/mf(mu_t, 0);
#define VOL_TRACE_LEAF \
  if(tmin > max_dist)\
  {\
    if(transmittance) *transmittance = mf_exp(mf_mul(mf_set1(-thickness), mu_t));\
    if(last_emission) *last_emission = mf_set1(0.f);\
    if(last_density)  *last_density = 0.f;\
    if(emission) *emission = mf_mul(*emission, *transmittance);\
    return FLT_MAX;\
  }\
  const float dist = _vol_trace_leaf(tree, mu_t, log_rand, rand, &thickness, result[0], result[1], tmin, MIN(max_dist,tmax), lambda, emission, transmittance, last_density, last_emission);\
  if(dist > max_dist) return FLT_MAX;\
  if(dist > 0.0f) return dist;
#define VOL_TRACE_MISS \
  if(transmittance) *transmittance = mf_set1(1.0f);\
  return FLT_MAX;
#define VOL_TRACE_RETURN \
  if(transmittance) *transmittance = mf_exp(mf_mul(mf_set1(-thickness), mu_t));\
  if(last_emission) *last_emission = mf_set1(0.f);\
  if(last_density) *last_density = 0.f;\
  if(emission) *emission = mf_mul(*emission, *transmittance);\
  return FLT_MAX;
#include "vol/trace_impl.inc"
}

// point lookup in world space
static inline int vol_lookup(vol_tree_t *tree, float i, float j, float k, vol_channel_t ch, vol_interpolation_t interpolation, int lod, float time, float *res)
{
  float pos[3] = {i, j, k};
  vol_transform_w2o(tree, pos, 1); // transform to object space
  if(interpolation & s_vol_smooth)
  {
    float size = tree->voxel_size * (1<<(3*lod));
    vol_interpolate_smooth(&pos[0], &pos[1], &pos[2], &time, size, 1.0/VOL_MOTION_SAMPLES);
  }
  pos[0] = (pos[0] - tree->aabb[0])/tree->voxel_size; // transform to voxel index space
  pos[1] = (pos[1] - tree->aabb[1])/tree->voxel_size;
  pos[2] = (pos[2] - tree->aabb[2])/tree->voxel_size;
  return vol_sample(tree, pos[0], pos[1], pos[2], lod, ch, time, res);
}

static inline float _vol_trace_leaf_line(
    const vol_tree_t *tree,
    const mf_t mu_t,
    const float log_rand,
    const float rand,
    float *thickness,
    float density,
    float temperature,
    const float tmin,
    const float tmax,
    const mf_t lambda,
    mf_t  *emission,
    mf_t  *transmittance,
    float *last_density,
    mf_t  *last_emission)
{
  float new_thickness = *thickness + density * MAX(0.0f, tmax-tmin);
  if(log_rand >= *thickness && log_rand < new_thickness)
  {
    float dist = 0;
    if(new_thickness > *thickness)
    {
      const float c0 = 1.0f-expf(-*thickness*mf(mu_t,0)), c1 = 1.0f-expf(-new_thickness*mf(mu_t,0));
      const float r2 = (rand - c0)/(c1-c0);
      const float M = (c1-c0)/(1.0-c0);
      dist = - logf(1.0f-r2*M)/(mf(mu_t,0)*density);
    }

    if (dist > 0 && dist < tmax-tmin)
    {
      const mf_t transmittance_voxel = mf_exp(mf_mul(mf_set1(-density * dist), mu_t));
      const mf_t Le = mf_mul(mf_set1(density), tree->shader(0, temperature, lambda));
      mf_t emission_voxel = density > 0.f ?
        // (*transmittance) * Le / (mu_t*density) * (1-transmittance_voxel)
        mf_mul(mf_div(mf_mul(*transmittance, Le), mf_mul(mu_t, mf_set1(density))), mf_sub(mf_set1(1.0f), transmittance_voxel))
        : 0.f;
      if (last_emission) *last_emission = Le;
      if (last_density)  *last_density = density;
      if (emission) *emission = mf_add(*emission, emission_voxel);
      *transmittance = mf_mul(*transmittance, transmittance_voxel);
      if (emission) *emission = mf_add(*emission,
          mf_div(mf_mul(Le, mf_sub(mf_set1(1.0f), *transmittance)), mf_mul(mu_t, mf_set1(density))));
      if (emission) *emission = mf_mul(*emission,
          mf_mul(mf_mul(mu_t, mf_set1(density)), *transmittance));
      return MAX(1e-6f, tmin + dist);
    }
  }

  const float dist = MAX(0.0f, tmax-tmin);
  const mf_t transmittance_voxel = mf_exp(mf_mul(mf_set1(-density * dist), mu_t));
  if(emission)
  {
    const mf_t Le = mf_mul(mf_set1(density), tree->shader(0, temperature, lambda));
    mf_t emission_voxel = density > 0.f ?
          mf_mul(*transmittance, mf_div(mf_mul(Le, mf_sub(mf_set1(1.0f), transmittance_voxel)), mf_mul(mu_t, mf_set1(density))))
      // (*transmittance) * Le / (mu_t*density) * (1-transmittance_voxel)
        : 0.f;
    *emission = mf_add(*emission, emission_voxel);
  }
  *thickness = new_thickness;
  *transmittance *= transmittance_voxel;
  return -1;
}

// pdf is transmittance (*mu_t at end vertex when ending in the volume)
static inline float vol_trace_sample_line(
    const vol_tree_t *tree,
    const float *pos_ws,
    const float *dir_ws,
    const float max_dist,             
    const mf_t mu_t,
    const float rand,                  
    vol_interpolation_t interpolation,  
    int lod,          
    float time,        
    const mf_t lambda, 
    mf_t  *emission,     
    mf_t  *transmittance, 
    mf_t  *last_emission,
    float *last_density) 
{
  assert(transmittance);
  float pos[3] = {pos_ws[0], pos_ws[1], pos_ws[2]};
  float dir[3] = {dir_ws[0], dir_ws[1], dir_ws[2]};
  vol_transform_w2o(tree, pos, 1);
  vol_transform_w2o(tree, dir, 0);
  // very slight speedup gained by optimising out a few exp by doing
  // cdf in log space (optical thickness).
  // M = 1-tau_voxel(tmax-tmin)
#define VOL_TRACE_INIT float thickness = 0.0f; *transmittance = mf_set1(1.f); if (emission) *emission = mf_set1(0.f);\
  const float log_rand = -logf(1.0f-rand)/mf(mu_t,0);
#define VOL_TRACE_LEAF \
  if(tmin >= max_dist)\
  {\
    if(transmittance) *transmittance = mf_exp(mf_mul(mf_set1(-thickness), mu_t));\
    if(last_emission) *last_emission = mf_set1(0.f);\
    if(last_density)  *last_density = 0.f;\
    if(emission) *emission = mf_mul(*emission, *transmittance);\
    return FLT_MAX;\
  }\
  const float dist = _vol_trace_leaf_line(tree, mu_t, log_rand, rand, &thickness,\
      result[0], result[1], tmin, MIN(max_dist, tmax), lambda, emission, transmittance, last_density, last_emission);\
  if(dist > 0.0 && dist < max_dist) return dist;
#define VOL_TRACE_MISS \
  return FLT_MAX;
#define VOL_TRACE_RETURN \
  if (last_emission) *last_emission = mf_set1(0.f);\
  if (last_density) *last_density = 0.f;\
  if (emission) *emission = mf_mul(*emission, *transmittance);\
  return FLT_MAX;
#include "vol/trace_impl.inc"
}
