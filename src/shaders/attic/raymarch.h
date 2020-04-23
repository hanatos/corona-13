#pragma once
// this doesn't perform woodcock tracking but works on the assumption that
// voxels represent piecewise homogeneous media. conceptually, this code
// constructs a cdf for every voxel along a ray, samples from that, and
// then samples a distance within this homogeneous block. in practice, the
// random number is chosen up front, so the cdf is never computed as whole.


#define MODE_T 0
#define MODE_SAMPLE 1
#define MODE_PDF 2
// MODE:
// T      - compute and return transmittance, use rdist as tmax
// SAMPLE - use rdist as random to compute scattering distance 
// PDF    - use rdist as distance to compute pdf
template<int MODE>
inline float volume_march(
    const medium_t *s,    // medium shader struct
    const path_t *p,      // light transport path
    const int e,          // edge number to march
    const float rdist)    // random number to sample the bin or distance
{
  float transmittance = 0.0f;
  float cdf = 0.0f;
  int i = 0;

  float start[3], end[3];
  const float tmin = volume_get_interval(p, e,
      MODE==MODE_T?rdist:FLT_MAX,
      s->ws_aabb_density, s->density->constTransform(), start, end);
  if(MODE==MODE_SAMPLE && tmin < 0.0f) return FLT_MAX;
  if(MODE==MODE_PDF    && tmin < 0.0f) return 1.0;
  if(MODE==MODE_T      && tmin < 0.0f) return 1.0;

  openvdb::FloatGrid::ConstAccessor a = s->density->getConstAccessor();
  openvdb::Vec3SGrid::ConstAccessor ma = s->velocity->getConstAccessor();
  const openvdb::math::Transform &at = s->density->constTransform();
  const openvdb::math::Transform &mat = s->velocity->constTransform();
  // step from start to end in one-voxel steps, using a simple dda.
  int dim = 0;
  if(fabsf(end[1] - start[1]) > fabsf(end[0] - start[0])) dim = 1;
  if(fabsf(end[2] - start[2]) > fabsf(end[dim] - start[dim])) dim = 2;
  float step[3];
  for(int k=0;k<3;k++) step[k] = end[k] - start[k];
  const float norm = fabsf(step[dim]);
  for(int k=0;k<3;k++)
  {
    step[k] /= norm;
    // if one of the steps collapses to nan, that means we divided by something silly small, i.e. start ~= end.
    if(MODE==MODE_SAMPLE && !(step[k] == step[k])) return FLT_MAX;  // fly right through
    if(MODE==MODE_PDF    && !(step[k] == step[k])) return 1.0;      // must have intersected geo, edge pdf is 1
    if(MODE==MODE_T      && !(step[k] == step[k])) return 1.0;      // transmittance == 1
  }
  // world space distance of one step:
  const float dist = s->scale * sqrtf(dotproduct(step,step)) * VOL_STEP_INDIRECT;
  int pdf_bin = (MODE==MODE_PDF) ? (rdist >= FLT_MAX ? INT_MAX : (rdist-tmin)/dist): 0.0;
  // FIXME: wtf stupid clang++ idiosyncracy :(
  if(MODE == MODE_PDF && pdf_bin < 0) pdf_bin = INT_MAX;
  assert(MODE!=MODE_PDF || pdf_bin >= 0);
  float curr[3] = {start[0], start[1], start[2]};
  // random offset to scramble ray marching bias:
  const float randlen = points_rand(rt.points, common_get_threadid());
  for(int k=0;k<3;k++) curr[k] += step[k]*randlen;
  int eyeray = 0;//(p->v[e-1].mode & s_sensor); // fails for light tracing :(
  while(1)
  {
    // quit loop:
    for(int k=0;k<3;k++)
    {
      if(step[k] > 0 && (int)curr[k] > s->ss_aabb_density[3+k]) goto break_outer;
      if(step[k] < 0 && (int)curr[k] < s->ss_aabb_density[k]) goto break_outer;
    }
    if(step[dim] > 0 && (int)curr[dim] > (int)end[dim]) goto break_outer;
    if(step[dim] < 0 && (int)curr[dim] < (int)end[dim]) goto break_outer;

    // switch off filtering after 50% of the light transmittance is gone:
    // if(eyeray && (transmittance > -logf(.5f)/dist)) eyeray = 0; // switched off this optimisation for fair comparison with woodcock tracking
    const float mu_t = s->sigma_t * volume_lookup_ss(curr, p->time, a, ma, at, mat, eyeray, s->motion_blur);

    transmittance += mu_t;

    if(MODE != MODE_T)
    {
      float new_cdf = cdf + (1.0-cdf)*(1.0-expf(-mu_t*dist));
      if(MODE == MODE_SAMPLE)
      {
        if((rdist >= cdf && rdist < new_cdf) || (new_cdf == 1.0))
        {
          assert(new_cdf-cdf > 0.0);
          const float r2 = (rdist - cdf)/(new_cdf-cdf); // rescale
          // TODO: use mu_t to do sub-voxel exponential distribution!
          assert(1e-6f + tmin + (i + r2)*dist > 0.0f);
          return 1e-6f + tmin + (i + r2)*dist;
        }
      }
      else if(MODE == MODE_PDF)
      {
        // TODO: should also multiply expf(-mu_t dist) instead of dividing out uniform distribution.
        //       the corresponding mu_t from normalisation can be found in form of rv in sample() and brdf().
        if(i == pdf_bin)
          return (new_cdf - cdf)/dist;
      }
      i++; // i belongs to cdf, not new_cdf
      cdf = new_cdf;
    }
    else // if MODE == MODE_T
    { // early out for transmittance < 1e-10
      if(transmittance > -logf(1e-10f)/dist) return expf(-transmittance*dist);
    }
    
    if(!eyeray)
      for(int k=0;k<3;k++) curr[k] += VOL_STEP_INDIRECT * step[k];
    else
      for(int k=0;k<3;k++) curr[k] += step[k];
  }
break_outer:
  // exhausted voxels while searching for cdf > rand
  if(MODE == MODE_SAMPLE)
    return FLT_MAX;
  // pdf to go through the whole volume is transmittance, too (MODE_T and MODE_PDF)
  return expf(-transmittance*dist);
}
