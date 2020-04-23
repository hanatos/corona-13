// this little tool converts openvdb .vdb files to hierarchical grid .vol files.
#include "vol/vol.h"
#include "vol/shaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

// if net set, it will use scatter (pointcloud style).
// scatter may have the potential to be faster (doesn't rasterise empty space)
// but gather turns out to have more consistent performance in practice
#define RAST_GATHER

int num_files = 0;
openvdb::FloatGrid::Ptr temperature;
openvdb::FloatGrid::Ptr density;
openvdb::Vec3SGrid::Ptr velocity;
float voxmul = 1; // 8 for lod
float max_vel[6];
float velocity_scale = 1.0f;

const int subsample_grid = 8; // subsampling ratio for coarse grid
float *velocity_grid = 0;
int    velocity_grid_size[3] = {0};
float  velocity_grid_aabb[6] = {0};

int get_max_velocity(const float *aabb, float *v)
{
  for(int k=0;k<6;k++) v[k] = 0.0f;
  if(!velocity_grid) return 0;
  // get int coords of corner points and then walk all of them:
  int m[3] = {0}, M[3] = {0};
  for(int k=0;k<3;k++)
  {
    m[k] = (int)CLAMP(velocity_grid_size[k] * (aabb[k]   - velocity_grid_aabb[k])/(velocity_grid_aabb[k+3]-velocity_grid_aabb[k]),       0, velocity_grid_size[k]-1);
    M[k] = (int)CLAMP(velocity_grid_size[k] * (aabb[k+3] - velocity_grid_aabb[k])/(velocity_grid_aabb[k+3]-velocity_grid_aabb[k]) + .5f, 0, velocity_grid_size[k]-1);
  }
  int empty = 1;
  for(int k=m[2];k<=M[2];k++)
  for(int j=m[1];j<=M[1];j++)
  for(int i=m[0];i<=M[0];i++)
  {
    const float *mv = velocity_grid + 7*(i + j*velocity_grid_size[0] + k*velocity_grid_size[0]*velocity_grid_size[1]);
    for(int k=0;k<3;k++) if(mv[k] < v[k]) v[k] = mv[k];
    for(int k=3;k<6;k++) if(mv[k] > v[k]) v[k] = mv[k];
    if(mv[6] != 0.0f) empty = 0;
  }
  return empty;
}

float xorshift_rng(uint64_t &state0, uint64_t &state1)
{
  uint64_t s1 = state0;
  uint64_t s0 = state1;
  state0 = s0;
  s1 ^= s1 << 23;
  s1 ^= s1 >> 17;
  s1 ^= s0;
  s1 ^= s0 >> 26;
  state1 = s1;
  return (state0 + state1) / ((double)((uint64_t)-1) + 1.0);
  // uint32_t v = 0x3f800000 | ((state0+state1)>>41); // faster than double version but breaks strict aliasing.
  // return (*(float*)&v) - 1.0f;
}

// callback used by vol interface to sample the tree
vol_payload_type_t payload_fill(
    void *data, vol_payload_uncompressed_t *payload,
    const float *aabb,
    int force_static)
{
#ifndef RAST_GATHER
  const int size[] = { 8,8,8 };
  // we want to have at least as many steps as output voxels along a velocity vector:
  const float velocity_step = MIN(MIN(
    (aabb[3] - aabb[0])/size[0],
    (aabb[4] - aabb[1])/size[1]),
    (aabb[5] - aabb[2])/size[2]);

  uint64_t state0 = aabb[0]+aabb[1]+aabb[2], state1 = 1;
  vol_payload_type_t have_data = s_vol_empty;

  // float samples[512][VOL_MOTION_SAMPLES] = {{0}};

  // grow input box by global motion bounds
  float input_aabb[6];

  for(int k=0;k<6;k++) input_aabb[k] = aabb[k] + max_vel[k];
  for(int k=0;k<2;k++)
  { // could iterate this a bit for even tighter bounds:
    float tighter_max_vel[6] = {0.0f};
    int empty = get_max_velocity(input_aabb, tighter_max_vel);
    if(empty) return s_vol_empty;
    for(int k=0;k<6;k++) input_aabb[k] = aabb[k] + tighter_max_vel[k];
  }

  // number of voxels to sample on interesting input block.  we clamp that to
  // avoid insane performance penalties for high motion blur, at the cost of a
  // blurry volume (pretend input voxels are actually bigger, and just
  // supersample during rasterisation)
  int input_size[3] = {
    // XXX this clamping thing resulted in grossly different results and performance :(
    // XXX if nothing else, it distorts the density of the volume, we need input/output voxel size compensation!
    // XXX at 30x the volume is nice and smooth! but 10x slower than a 5x multiplier
    CLAMP((int)(.5f + size[0]*(input_aabb[3]-input_aabb[0])/(aabb[3]-aabb[0])), size[0], 8*size[0]),
    CLAMP((int)(.5f + size[1]*(input_aabb[4]-input_aabb[1])/(aabb[4]-aabb[1])), size[1], 8*size[1]),
    CLAMP((int)(.5f + size[2]*(input_aabb[5]-input_aabb[2])/(aabb[5]-aabb[2])), size[2], 8*size[2])
  };
  const float voxsize[3] = { // world space voxel size in block we consider on input
    (input_aabb[0+3]-input_aabb[0])/input_size[0],
    (input_aabb[1+3]-input_aabb[1])/input_size[1],
    (input_aabb[2+3]-input_aabb[2])/input_size[2]
  };
  const float voxsize_out[3] = { // world space voxel size of output voxels
    (aabb[0+3]-aabb[0])/size[0],
    (aabb[1+3]-aabb[1])/size[1],
    (aabb[2+3]-aabb[2])/size[2]
  };
  const int supersample = // adjust grid sampling densities
    // note that vdb will be read out only once per voxel.
    CLAMP(
        voxsize[0]/voxsize_out[0]*
        voxsize[1]/voxsize_out[1]*
        voxsize[2]/voxsize_out[2],
    // XXX this is wrong. it should be wrt voxel sizes, not voxel counts (due to clamping above these may be different things)
        //input_size[0]/(float)size[0]*
        //input_size[1]/(float)size[1]*
        //input_size[2]/(float)size[2],
        1, 32);

  // if(input_size[0] > 8 || input_size[1] > 8 || input_size[2] > 8)
    // fprintf(stderr, "input size %d %d %d\n", input_size[0], input_size[1], input_size[2]);
  for(int ik=0;ik<input_size[2];ik++)
  for(int ij=0;ij<input_size[1];ij++)
  for(int ii=0;ii<input_size[0];ii++)
  {
    // sample input, fill leaf struct
    float x[3];
    x[0] = input_aabb[0] + voxsize[0]*(ii+.5f);
    x[1] = input_aabb[1] + voxsize[1]*(ij+.5f);
    x[2] = input_aabb[2] + voxsize[2]*(ik+.5f);
    openvdb::Vec3d world(x[0], x[1], x[2]);
    // get payload at this index:
    const float dens = density->tree().getValue(density->constTransform().worldToIndexCellCentered(world));
    if(dens <= 0.0) continue;
    float temp = 0.0f;
    openvdb::Vec3s v(0,0,0);
    if(temperature)
      temp = temperature->tree().getValue(temperature->constTransform().worldToIndexCellCentered(world));
    if(velocity)
      v = velocity_scale * velocity->tree().getValue(velocity->constTransform().worldToIndexCellCentered(world));

    const float vl = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    float tmin = 0;   // time interval of intersection
    float tmax = 1.0f;
#if 1
    for(int k=0;k<3;k++)
    {
      // grow box a bit, we're only testing the centre point
      const float t0 = (aabb[k+3] + .5f * voxsize[k] - x[k])/v[k];
      const float t1 = (aabb[k]   - .5f * voxsize[k] - x[k])/v[k];
      if(t0 <= t1)
      {
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
      }
      else if(t1 < t0)
      {
        tmin = t1 > tmin ? t1 : tmin;
        tmax = t0 < tmax ? t0 : tmax;
      }
    }
    // miss bounding box?
    if(tmin > tmax) continue;
#endif

    // splat a number of sub-input voxel samples.
    for(int r=0;r<supersample;r++)
    {
      float pos0[3], pos1[3];
      pos0[0] = input_aabb[0] + (ii+xorshift_rng(state0, state1))/input_size[0] * (input_aabb[3+0]-input_aabb[0]);
      pos0[1] = input_aabb[1] + (ij+xorshift_rng(state0, state1))/input_size[1] * (input_aabb[3+1]-input_aabb[1]);
      pos0[2] = input_aabb[2] + (ik+xorshift_rng(state0, state1))/input_size[2] * (input_aabb[3+2]-input_aabb[2]);
      for(int k=0;k<3;k++) pos1[k] = pos0[k] + v[k];

      const int vsamples = MAX(1.0f, (tmax-tmin)*vl / velocity_step);
      // round time to nearest bins outwards:
      const float min_time = ((int)(tmin * VOL_MOTION_SAMPLES))/(float)VOL_MOTION_SAMPLES;
      const float max_time = ((int)(tmax * VOL_MOTION_SAMPLES + 0.99))/(float)VOL_MOTION_SAMPLES;
      // fprintf(stderr, "sampling %g %g with %d\n", min_time, max_time, vsamples);
      for(int s=0;s<vsamples;s++)
      {
        const float frac = s/(float)vsamples;
        const float next_frac = (s+1.0f)/(float)vsamples;
        const float time = MAX(0.0f, min_time * (1.0f-frac) + max_time * frac);
        const float next_time = MIN(1.0f, min_time * (1.0f-next_frac) + max_time * next_frac);

        float pos[3];
        for(int k=0;k<3;k++) pos[k] = (1.0f-time) * pos0[k] + time * pos1[k];
        vol_index_t vxo = {0};
        const int i = (int)((pos[0] - aabb[0])/(aabb[3]-aabb[0]) * size[0]);
        const int j = (int)((pos[1] - aabb[1])/(aabb[4]-aabb[1]) * size[1]);
        const int k = (int)((pos[2] - aabb[2])/(aabb[5]-aabb[2]) * size[2]);
        if(i < 0 || i >= size[0] ||
           j < 0 || j >= size[1] ||
           k < 0 || k >= size[2]) continue;
        vxo.i = i; vxo.j = j; vxo.k = k;

        // const int t = (int)CLAMP(time * VOL_MOTION_SAMPLES, 0, VOL_MOTION_SAMPLES-1);

        for(int t=time * VOL_MOTION_SAMPLES; t<next_time * VOL_MOTION_SAMPLES;t++)
        {
          // write:
          payload->d[vxo.idx][t] = dens;
          payload->t[vxo.idx][t] = temp;
          // average result:
          // const float wo = samples[vxo.idx][t];
          // payload->d[vxo.idx][t] = wo/(wo+1.0f) * payload->d[vxo.idx][t] + 1.0f/(wo+1.0f) * dens;
          // payload->t[vxo.idx][t] = wo/(wo+1.0f) * payload->t[vxo.idx][t] + 1.0f/(wo+1.0f) * temp;
          // samples[vxo.idx][t]++;
          // just normalise to allow changing density over time! (contract/expand)
          // we don't have enough samples, not random enough. leads to artifacts:
          // payload->d[vxo.idx][t] += 1.0f/supersample * dens;
          // payload->t[vxo.idx][t] += 1.0f/supersample * temp;
        }
        have_data = force_static ? s_vol_static : s_vol_full;
      }
    }
  }
  return have_data;

#else
  vol_payload_type_t have_data = s_vol_empty;

  // grow input box by global motion bounds, early out:
  float input_aabb[6];
  for(int k=0;k<6;k++) input_aabb[k] = aabb[k] + max_vel[k];
  for(int k=0;k<2;k++)
  { // could iterate this a bit for even tighter bounds:
    float tighter_max_vel[6] = {0.0f};
    int empty = get_max_velocity(input_aabb, tighter_max_vel);
    if(empty) return s_vol_empty;
    for(int k=0;k<6;k++) input_aabb[k] = aabb[k] + tighter_max_vel[k];
  }

  // didn't do much difference in performance at all:
  // openvdb::Vec3SGrid::ConstAccessor accv = velocity->getConstAccessor();
  // openvdb::FloatGrid::ConstAccessor accd = density->getConstAccessor();
  // openvdb::FloatGrid::ConstAccessor acct = temperature->getConstAccessor();

  for(int i=0;i<512;i++)
  {
    vol_index_t vx = {0};
    vx.idx = i;
    // sample input, fill leaf struct
    float x[3];
    x[0] = aabb[0] + (aabb[0+3]-aabb[0])*(vx.i+.5f)/8.0f;
    x[1] = aabb[1] + (aabb[1+3]-aabb[1])*(vx.j+.5f)/8.0f;
    x[2] = aabb[2] + (aabb[2+3]-aabb[2])*(vx.k+.5f)/8.0f;
    openvdb::Vec3d world(x[0], x[1], x[2]);
#if 1 // motion blur
    if(velocity)
    {
      openvdb::math::Coord vpos = velocity->constTransform().worldToIndexCellCentered(world);
      // openvdb::Vec3f v = - velocity_scale * accv.getValue(vpos);
      openvdb::Vec3f v = - velocity_scale * velocity->tree().getValue(vpos);
      openvdb::Vec3f world_end = world + v;
      const float vlen = sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/velocity_scale; // length of motion vector in voxels
      const int vsamples = MIN(VOL_MOTION_SAMPLES, (int)vlen + 1);
      for(int s=0;s<vsamples;s++)
      {
        float time = s/(vsamples-1.0f);
        // resample to num_files and fractional time therein:
        // float ftime = num_files * time;
        // int f = CLAMP(ftime, 0, num_files-1);
        float fractime = time ;// ftime - f;
        openvdb::Vec3f advected = (1-fractime) * world + fractime * world_end;
        float temp = 0.0f;
        if(temperature)
          temp = temperature->tree().getValue(temperature->constTransform().worldToIndexCellCentered(advected));
          // temp = acct.getValue(density->constTransform().worldToIndexCellCentered(advected));
        const float dens = density->tree().getValue(density->constTransform().worldToIndexCellCentered(advected));
        // const float dens = accd.getValue(density->constTransform().worldToIndexCellCentered(advected));
        if(dens > 0.0) have_data = force_static ? s_vol_static : s_vol_full;
        for(int t=s/(float)vsamples * VOL_MOTION_SAMPLES;
                t<(s+1)/(float)vsamples * VOL_MOTION_SAMPLES;t++)
        {
          assert(t >= 0);
          assert(t < VOL_MOTION_SAMPLES);
          payload->d[i][t] = dens;
          payload->t[i][t] = temp;
        }
      }
    }
    else
#endif // no motion blur, to debug compression
    {
      const float dens = density->tree().getValue(density->constTransform().worldToIndexCellCentered(world));
      const float temp = temperature->tree().getValue(temperature->constTransform().worldToIndexCellCentered(world));
      for(int t=0;t<VOL_MOTION_SAMPLES;t++)
      {
        payload->d[i][t] = dens;
        payload->t[i][t] = temp;
        if(payload->d[i][t] > 0.0) have_data = s_vol_static;
      }
    }
  }
  return have_data;
#endif
}

int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "[ovdb2vol] usage: %s input[%%04d].vdb [-s, --static] [--shader-id x] [lod]\n", arg[0]);
    exit(1);
  }
  openvdb::initialize(); // multiple calls are okay

  int animate = 0;
  // XXX disabled with new rasterisation
  // if(strstr(arg[1], "%04d")) animate = 1;
  char filename[2048];
  int beg = 9999, end = 9999;

  int motion_blur = 1;
  vol_shader_id_t shader_id = s_vol_shader_blackbody;
  if(argc > 2) for(int a=2;a<argc;a++)
  {
    if(arg[a][0] == '-' && (arg[a][1] == 's' || !strcmp(arg[a], "--static")))
      motion_blur = 0;
    else if(arg[a][0] == '-' && !strcmp(arg[a], "--shader-id") && (++a < argc))
      shader_id = vol_shader_id_t(atol(arg[a]));
    else voxmul = powf(8.0, atol(arg[a]));
  }

  // find %04d in file name and if so loop the following until file not found:
  num_files = 0;
  if(animate) for(int i=0;i<=9999;i++)
  {
    snprintf(filename, sizeof(filename), arg[1], i);
    struct stat sb;
    memset(&sb, 0, sizeof(sb));
    stat(filename, &sb);
    if(S_ISREG(sb.st_mode))
    {
      if(beg == 9999) beg = i;
      end = i;
    }
    else if(beg < 9999) break;
  }
  if(animate)
    fprintf(stderr, "[ovdb2vol] processing %d frames %d--%d\n", end-beg+1, beg, end);

  // temperature = new openvdb::FloatGrid::Ptr[end-beg+1];
  // density     = new openvdb::FloatGrid::Ptr[end-beg+1];
  // velocity    = new openvdb::Vec3SGrid::Ptr[end-beg+1];
  
  const double begin_time = common_time_wallclock();
  for(int i=beg;i<=end;i++)
  {
    if(animate) snprintf(filename, sizeof(filename), arg[1], i);
    else        strncpy(filename, arg[1], sizeof(filename));

    openvdb::io::File file(filename);
    try
    {
      file.open(); // read header only
    }
    catch (openvdb::IoError &e)
    {
      fprintf(stderr, "[ovdb2vol] could not open file `%s'!\n", filename);
      return 1;
    }

    for(openvdb::io::File::NameIterator it = file.beginName(); it != file.endName(); ++it)
    {
      if(it.gridName() == "density")
        // density[num_files] = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(it.gridName()));
        density = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(it.gridName()));
      if(it.gridName() == "temperature" || it.gridName() == "emission")
        temperature = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(it.gridName()));
      if(motion_blur &&
        (it.gridName() == "v" || it.gridName() == "velocity" || it.gridName() == "vel"))
        velocity = openvdb::gridPtrCast<openvdb::Vec3SGrid>(file.readGrid(it.gridName()));
    }
    file.close();
    num_files++;
  }

  // assume last bounding box is most expanded (fireballs/explosions)
  openvdb::CoordBBox bbox = density->evalActiveVoxelBoundingBox();
  openvdb::Coord dim = density->evalActiveVoxelDim();
  std::cout << "[ovdb2vol] grid dimensions: " << dim << std::endl;

  openvdb::Vec3d wm = density->constTransform().indexToWorld(bbox.min());
  openvdb::Vec3d wM = density->constTransform().indexToWorld(bbox.max());
  const float worldspace[6] = {(float)wm[0], (float)wm[1], (float)wm[2], (float)wM[0], (float)wM[1], (float)wM[2]}; 
  const float voxel_size = voxmul * fminf((wM[0]-wm[0])/dim[0], fminf((wM[1]-wm[1])/dim[1], (wM[2]-wm[2])/dim[2]));

  // velocity is given as voxel offset, we need it in world space:
  velocity_scale = (wM[0] - wm[0]) / dim[0];
  if(!motion_blur) velocity_scale = 0.0f;

  // init max velocity used for coarse motion bound
  for(int k=0;k<6;k++) max_vel[k] = 0.0f;
  velocity_grid = 0;
  if(velocity)
  {
    openvdb::Vec3SGrid::ConstAccessor a = velocity->getConstAccessor();
    openvdb::CoordBBox bboxv = velocity->evalActiveVoxelBoundingBox();
    openvdb::Coord dimv = velocity->evalActiveVoxelDim();
    // create velocity grid:
    openvdb::Vec3d vwm = velocity->constTransform().indexToWorld(bboxv.min());
    openvdb::Vec3d vwM = velocity->constTransform().indexToWorld(bboxv.max());
    velocity_grid_aabb[0] = vwm[0];
    velocity_grid_aabb[1] = vwm[1];
    velocity_grid_aabb[2] = vwm[2];
    velocity_grid_aabb[3] = vwM[0];
    velocity_grid_aabb[4] = vwM[1];
    velocity_grid_aabb[5] = vwM[2];
    velocity_grid_size[0] = (dimv[0]-1)/subsample_grid+1;
    velocity_grid_size[1] = (dimv[1]-1)/subsample_grid+1;
    velocity_grid_size[2] = (dimv[2]-1)/subsample_grid+1;
    velocity_grid = (float *)malloc(sizeof(float)*7*velocity_grid_size[0]*velocity_grid_size[1]*velocity_grid_size[2]);
    std::cout << "[ovdb2vol] creating an auxiliary "
      <<velocity_grid_size[0]<<"x" <<velocity_grid_size[1]<<"x" <<velocity_grid_size[2]<<" velocity grid.." << std::endl;
    memset(velocity_grid, 0, sizeof(float)*7*velocity_grid_size[0]*velocity_grid_size[1]*velocity_grid_size[2]);
#pragma omp parallel for schedule(static) collapse(3) firstprivate(a) default(shared)
    for(int k=0;k<dimv[2];k++) for(int j=0;j<dimv[1];j++) for(int i=0;i<dimv[0];i++)
    {
      openvdb::Vec3d world = velocity->constTransform().indexToWorld(openvdb::Coord(i, j, k));
      float d = density->tree().getValue(density->constTransform().worldToIndexCellCentered(world));
      if(d > 0.0)
      { // only count velocity if it actually moves any mass
        openvdb::Vec3s v = - velocity_scale * a.getValue(openvdb::Coord(i, j, k));
        for(int k=0;k<3;k++) if(v[k] < max_vel[k])   max_vel[k]   = v[k];
        for(int k=0;k<3;k++) if(v[k] > max_vel[3+k]) max_vel[3+k] = v[k];
        // get from worldspace coords just same as above when reading it:
        int ii = (int)CLAMP(velocity_grid_size[0] * (world[0]   - velocity_grid_aabb[0])/(velocity_grid_aabb[0+3]-velocity_grid_aabb[0]),       0, velocity_grid_size[0]-1);
        int jj = (int)CLAMP(velocity_grid_size[1] * (world[1]   - velocity_grid_aabb[1])/(velocity_grid_aabb[1+3]-velocity_grid_aabb[1]),       0, velocity_grid_size[1]-1);
        int kk = (int)CLAMP(velocity_grid_size[2] * (world[2]   - velocity_grid_aabb[2])/(velocity_grid_aabb[2+3]-velocity_grid_aabb[2]),       0, velocity_grid_size[2]-1);
        float *mv = velocity_grid + 7*(ii+velocity_grid_size[0]*jj + velocity_grid_size[0]*velocity_grid_size[1]*kk);
        for(int k=0;k<3;k++) if(v[k] < mv[k])   mv[k]   = v[k];
        for(int k=0;k<3;k++) if(v[k] > mv[3+k]) mv[3+k] = v[k];
        mv[6] = 1.0f;
      }
    }
  }
  std::cout << "[ovdb2vol] motion bounds: " << max_vel[0] << " " << max_vel[1] << " " << max_vel[2] <<
                                        " " << max_vel[3] << " " << max_vel[4] << " " << max_vel[5] << std::endl;
  const double mid_time = common_time_wallclock();

  char output[512];
  if(animate) snprintf(output, 512, arg[1], 0);
  else        snprintf(output, 512, "%s", arg[1]);
  char *c = output + strlen(output);
  for(;c > output && *c != '.';c--);
  snprintf(c, 5, ".vol");

  if(vol_create_tree(output, worldspace, voxel_size, 0, 0,
        &payload_fill, 0, shader_id, !motion_blur))
    exit(1);
  const double end_time = common_time_wallclock();
  std::cout << "[ovdb2vol] ingested vdb in "<<(mid_time-begin_time)<<"s, created tree in "<<(end_time-mid_time)<<"s."<<std::endl;
  free(velocity_grid);
  return 0;
}
