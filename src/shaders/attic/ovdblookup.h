#pragma once

inline void volume_get_interval_bounds(
    const float *pos,
    const float *dir,
    const float *aabb,
    float *tmin,
    float *tmax)
{
  // intersect aabb
  for(int k=0;k<3;k++)
  {
    const float t0 = (aabb[k+3] - pos[k])/dir[k];
    const float t1 = (aabb[k]   - pos[k])/dir[k];
    if(t0 <= t1)
    {
      *tmin = t0 > *tmin ? t0 : *tmin;
      *tmax = t1 < *tmax ? t1 : *tmax;
    }
    else
    {
      *tmin = t1 > *tmin ? t1 : *tmin;
      *tmax = t0 < *tmax ? t0 : *tmax;
    }
  }
}

// returns the start/end interval in (fractional) voxel index space coordinates,
// given the path edge p->e[e] and the world space voxel bounding box aabb
inline float volume_get_interval(
    const path_t *p,
    const int e,
    float tmax,
    const float *aabb,
    const openvdb::math::Transform &tf,
    float *start,
    float *end)
{
  // intersect aabb
  float tmin = 0.f;
  volume_get_interval_bounds(p->v[e-1].hit.x, p->e[e].omega, aabb, &tmin, &tmax);
  // can only work if start != end because only then we get a valid step direction:
  if(tmin < tmax)
  {
    // compute world space start/end point. already subtract aabb min:
    for(int k=0;k<3;k++) start[k] = p->v[e-1].hit.x[k] + p->e[e].omega[k] * tmin;
    for(int k=0;k<3;k++) end  [k] = p->v[e-1].hit.x[k] + p->e[e].omega[k] * tmax;
    // convert to sample space (in voxel grid):
    openvdb::Vec3d tmps = tf.worldToIndex(openvdb::Vec3d(start[0], start[1], start[2]));
    openvdb::Vec3d tmpe = tf.worldToIndex(openvdb::Vec3d(end[0], end[1], end[2]));
    for(int k=0;k<3;k++) start[k] = tmps[k];
    for(int k=0;k<3;k++) end[k] = tmpe[k];
    assert(tmin >= 0);
    return tmin;
  }
  return -1;
}

inline float volume_lookup_ss(
    const float *pos,                          // lookup coordinates in local index space of accessor a
    const float time,                          // path's time for motion blur
    openvdb::FloatGrid::ConstAccessor &a,      // value to return
    openvdb::Vec3SGrid::ConstAccessor &ma,     // motion vector field
    const openvdb::math::Transform &ta,        // transform of value grid
    const openvdb::math::Transform &tma,       // transform of motion vector grid
    const int eyeray,                          // flag to signify high quality interpolation or not
    const int motion_blur)                     // indicates we have a velocity field
{
  if(!motion_blur)
  { // no motion blur
    // triquadratic:
    if(eyeray)
    {
      float value;
      // openvdb::tools::QuadraticSampler::sample(a, openvdb::Vec3d(pos[0], pos[1], pos[2]), value);
      openvdb::tools::BoxSampler::sample(a, openvdb::Vec3d(pos[0], pos[1], pos[2]), value);
      if(value < VOL_MU_T_CUTOFF) return 0.0f;
      return fmaxf(0.0, value); // interpolation might overshoot
    }
    else // unfiltered, but cached:
    {
      const float value = a.getValue(openvdb::Coord(pos[0], pos[1], pos[2]));
      if(value < VOL_MU_T_CUTOFF) return 0.0f;
      return value;
    }
  }
  else
  { // with motion blur:
    if(eyeray)
    {
      // triquadratic:
      openvdb::Vec3s v;
      openvdb::Vec3d world = ta.indexToWorld(openvdb::Vec3d(pos[0], pos[1], pos[2]));
      openvdb::tools::QuadraticSampler::sample(ma, tma.worldToIndex(world), v);
      openvdb::Vec3d advected;
      for(int k=0;k<3;k++) advected[k] = pos[k] + v[k]*(time - .5f);
      float value;
      // openvdb::tools::QuadraticSampler::sample(a, advected, value);
      // return fmaxf(0.0, value); // interpolation might overshoot
      openvdb::tools::BoxSampler::sample(a, advected, value);
      if(value < VOL_MU_T_CUTOFF) return 0.0f;
      return value;
    }
    else
    { // unfiltered:
      openvdb::Vec3d world = ta.indexToWorld(openvdb::Vec3d(pos[0], pos[1], pos[2]));
      openvdb::Vec3s v = ma.getValue(openvdb::Coord(tma.worldToIndexCellCentered(world)));

      // fprintf(stderr, "veloc %g %g %g * %g\n", v[0], v[1], v[2], time);
      float advected[3];
      for(int k=0;k<3;k++) advected[k] = pos[k] + v[k]*(time - .5);
      const float value = a.getValue(openvdb::Coord(advected[0], advected[1], advected[2]));
      if(value < VOL_MU_T_CUTOFF) return 0.0f;
      return value;
    }
  }
}

inline float volume_lookup_ws(
    const float *pos,                          // lookup coordinates in world space
    const float time,                          // path's time for motion blur
    openvdb::FloatGrid::ConstAccessor &a,      // value to return
    openvdb::Vec3SGrid::ConstAccessor &ma,     // motion vector field
    const openvdb::math::Transform &ta,        // transform of value grid
    const openvdb::math::Transform &tma,       // transform of motion vector grid
    const int eyeray,                          // flag to signify high quality interpolation or not
    const int motion_blur)
{
  const openvdb::Vec3d v = ta.worldToIndex(openvdb::Vec3d(pos[0], pos[1], pos[2]));
  float sspos[3] = {(float)v[0], (float)v[1], (float)v[2]};
  return volume_lookup_ss(sspos, time, a, ma, ta, tma, eyeray, motion_blur);
}

