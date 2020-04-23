#pragma once

#include "points.h"

// called `delta tracking' in [Novak et al. 2014].
static inline float woodcock_track(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float mu_t_max,    // majorant mu_t
    const float min_dist,    // starting point (distance)
    const float max_dist,    // max distance to track to
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  openvdb::FloatGrid::ConstAccessor a = s->density->getConstAccessor();
  openvdb::Vec3SGrid::ConstAccessor ma = s->velocity->getConstAccessor();
  const openvdb::math::Transform &at = s->density->constTransform();
  const openvdb::math::Transform &mat = s->velocity->constTransform();

  float curr[3];
  float dist = min_dist;
  while(1)
  {
    const float x1 = points_rand(rt.points, common_get_threadid());
    // XXX this log is apparently 16% of total runtime:
    dist -= logf(1.0f - x1)/mu_t_max;
    if(dist > max_dist) return FLT_MAX;
    const float x2 = points_rand(rt.points, common_get_threadid());
    for(int k=0;k<3;k++) curr[k] = pos[k] + dist*dir[k];
    // XXX this lookup is another 20% if i read this right
    const float mu_t = s->sigma_t * volume_lookup_ws(curr, time, a, ma, at, mat, eyeray, s->motion_blur); // XXX replace by ss lookup
    if(x2 <= mu_t/mu_t_max) break;
  }
  return dist;
}

// stupid way of computing transmittance with woodcock tracking,
// only one estimate and return 1 or 0.
static inline float woodcock_track_transmittance(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float mu_t_max,    // majorant mu_t
    const float max_dist,    // max distance to track to
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  const float dist = woodcock_track(s, pos, dir, mu_t_max, 0.0f, max_dist, time, eyeray);
  return dist > max_dist ? 1.0 : 0.0;
}

// this is `residual ratio tracking' [Novak et al. 2014].
static inline float woodcock_transmittance(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float mu_t_ctr,    // control variate mu_t (average)
    const float mu_t_rmx,    // residual majorant mu_t
    const float min_dist,    // start distance
    const float max_dist,    // end distance to evaluate
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  openvdb::FloatGrid::ConstAccessor a = s->density->getConstAccessor();
  openvdb::Vec3SGrid::ConstAccessor ma = s->velocity->getConstAccessor();
  const openvdb::math::Transform &at = s->density->constTransform();
  const openvdb::math::Transform &mat = s->velocity->constTransform();

  float curr[3] = {pos[0], pos[1], pos[2]};
  float d = min_dist;
  const float t_ctr = expf(-mu_t_ctr * (max_dist-min_dist));
  // if(t_ctr < VOL_MIN_TRANSMITTANCE) return 0.0f;
  float t_res = 1.0f;
  while(1)
  {
    const float x1 = points_rand(rt.points, common_get_threadid());
    d -= logf(1.0f - x1)/mu_t_rmx;
    if(d > max_dist) break;

    for(int k=0;k<3;k++) curr[k] = pos[k] + d*dir[k];
    const float mu_t = s->sigma_t * volume_lookup_ws(curr, time, a, ma, at, mat, eyeray, s->motion_blur); // XXX replace by ss lookup
    // update residual transmittance
    t_res *= 1.0f - (mu_t - mu_t_ctr)/mu_t_rmx;
    // if(t_res < VOL_MIN_TRANSMITTANCE) return 0.0f;
  }
  return t_ctr * t_res;
}


// return the pdf of a ray reaching a certain distance.
static inline float woodcock_track_pdf(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float mu_t_max,    // majorant mu_t
    const float max_dist,    // max distance to track to
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  const float dist = woodcock_track(s, pos, dir, mu_t_max, 0.0f, max_dist, time, eyeray);
  return dist > max_dist;
}

static inline float woodcock_track_grid(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float max_dist,    // max distance to track to
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  const float *aabb = s->ws_aabb_density;
  int p[3];
  float t[3], invdir[3];
  for(int k=0;k<3;k++) invdir[k] = 1.0f/dir[k];
  const int step[3] = {invdir[0] > 0 ? 1 : -1, invdir[1] > 0 ? 1 : -1, invdir[2] > 0 ? 1 : -1 };

  // intersect bounding box, get entry and exit distances tmin, tmax
  float tmin = 0.0f;
  float tmax = max_dist;
  volume_get_interval_bounds(pos, dir, aabb, &tmin, &tmax);
  // miss bounding box?
  if(tmin > tmax) return FLT_MAX;

  // size of voxel in world space:
  const float cellsize[3] = {(aabb[3+0]-aabb[0])/s->grid_size[0],(aabb[3+1]-aabb[1])/s->grid_size[1],(aabb[3+2]-aabb[2])/s->grid_size[2]};
  const float supervoxelwd = 1<<s->grid_shift;

  // entry point as grid index p[]:
  for(int k=0;k<3;k++) t[k] = pos[k] + tmin*dir[k]; // world space entry point
  for(int k=0;k<3;k++) p[k] = (int)fmaxf(0.0f, fminf(s->grid_size[k]-1e-5f, (t[k] - aabb[k])/cellsize[k]));
  for(int k=0;k<3;k++) p[k] >>= s->grid_shift; // now in supervoxel coordinates

  // init t[] to distance offsets to next voxel in every dimension by intersecting a voxel-sized aabb:
  for(int k=0;k<3;k++) 
    t[k] = (aabb[k] + (p[k]+((step[k]==1)&1)) * supervoxelwd * cellsize[k] - pos[k])*invdir[k];

  // step delta for t[] in each dimension:
  const float delta[3] = {cellsize[0]*supervoxelwd*fabsf(invdir[0]), cellsize[1]*supervoxelwd*fabsf(invdir[1]), cellsize[2]*supervoxelwd*fabsf(invdir[2]) };

  // int counter = 0;
  while(1)
  {
    // compute entry point distance t[mind] to next voxel:
    int mind = 0;
    if(t[2] < t[1] && t[2] < t[0]) mind = 2;
    else if(t[1] < t[0]) mind = 1;

    const float mu_max = s->grid[p[0] + s->grid_coarse_size[0] * (p[1] + s->grid_coarse_size[1]*p[2])].mu_max;
    if(mu_max > 0.0)
    {
      const float dist = woodcock_track(s, pos, dir, mu_max, tmin, t[mind], time, eyeray);
      if(dist < t[mind]) return dist;
    }

    tmin = t[mind];

    t[mind] += delta[mind];
    p[mind] += step[mind];
    if(p[mind] >= s->grid_coarse_size[mind] || p[mind] < 0) return FLT_MAX;
  }
}

static inline float woodcock_transmittance_grid(
    const medium_t *s,       // struct holding all vdb related data
    const float pos[3],      // origin of ray
    const float dir[3],      // direction of ray
    const float max_dist,    // max distance to track to
    const float time,        // time (motion blur)
    int eyeray)              // eye ray or not (interpolate voxels?)
{
  const float *aabb = s->ws_aabb_density;
  int p[3];
  float t[3], invdir[3];
  for(int k=0;k<3;k++) invdir[k] = 1.0f/dir[k];
  const int step[3] = {invdir[0] > 0 ? 1 : -1, invdir[1] > 0 ? 1 : -1, invdir[2] > 0 ? 1 : -1 };

  // intersect bounding box, get entry and exit distances tmin, tmax
  float tmin = 0.0f;
  float tmax = max_dist;
  volume_get_interval_bounds(pos, dir, aabb, &tmin, &tmax);
  // miss bounding box?
  if(tmin > tmax) return 1.0f;

  // size of voxel in world space:
  const float cellsize[3] = {(aabb[3+0]-aabb[0])/s->grid_size[0],(aabb[3+1]-aabb[1])/s->grid_size[1],(aabb[3+2]-aabb[2])/s->grid_size[2]};
  const float supervoxelwd = 1<<s->grid_shift;

  // entry point as grid index p[]:
  for(int k=0;k<3;k++) t[k] = pos[k] + tmin*dir[k]; // world space entry point
  for(int k=0;k<3;k++) p[k] = (int)fmaxf(0.0f, fminf(s->grid_size[k]-1e-5f, (t[k] - aabb[k])/cellsize[k]));
  for(int k=0;k<3;k++) p[k] >>= s->grid_shift; // now in supervoxel coordinates

  // init t[] to distance offsets to next voxel in every dimension by intersecting a voxel-sized aabb:
  for(int k=0;k<3;k++) 
    t[k] = (aabb[k] + (p[k]+((step[k]==1)&1)) * supervoxelwd * cellsize[k] - pos[k])*invdir[k];

  // step delta for t[] in each dimension:
  const float delta[3] = {cellsize[0]*supervoxelwd*fabsf(invdir[0]), cellsize[1]*supervoxelwd*fabsf(invdir[1]), cellsize[2]*supervoxelwd*fabsf(invdir[2]) };

  float transmittance = 1.0f;
  while(1)
  {
    // compute entry point distance t to next voxel:
    int mind = 0;
    if(t[2] < t[1] && t[2] < t[0]) mind = 2;
    else if(t[1] < t[0]) mind = 1;

    const float mu_c = s->grid[p[0] + s->grid_coarse_size[0] * (p[1] + s->grid_coarse_size[1]*p[2])].mu_c;
    if(mu_c > 0.0)
    {
      const float mu_r = s->grid[p[0] + s->grid_coarse_size[0] * (p[1] + s->grid_coarse_size[1]*p[2])].mu_r;
      transmittance *= woodcock_transmittance(s, pos, dir, mu_c, mu_r, tmin, t[mind], time, eyeray);
      if(transmittance < VOL_MIN_TRANSMITTANCE) return 0.0f;
    }

    tmin = t[mind];
    // if(t[mind] > tmax) return transmittance; one check is enough. need to benchmark which is faster
    t[mind] += delta[mind];
    p[mind] += step[mind];
    if(p[mind] >= s->grid_coarse_size[mind] || p[mind] < 0) return transmittance;
  }
}

