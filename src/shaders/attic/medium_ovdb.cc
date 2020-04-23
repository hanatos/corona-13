/*
    This file is part of corona-13.

    copyright (c) 2015 johannes hanika.

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

// heterogeneous medium, wrapping around openvdb.
// 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>


// need to compile c++ to link openvdb, but corona opens our callbacks as c11.
// clang whines about templates with c linkage, so we individually define the callbacks as extern "C".

#ifdef __cplusplus
extern "C" {
#endif
#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum.h"
#include "prims.h"
#ifdef __cplusplus
}
#endif

// volume positions with lower mu_t concentration are considered empty:
#define VOL_MU_T_CUTOFF 1e-5f
#define VOL_MIN_TRANSMITTANCE 1e-5f

#define VOL_METHOD 1 // ray marching
// #define VOL_METHOD 2 // woodcock tracking/residual ratio transmittance
// #define VOL_METHOD 3 // woodcock tracking/residual ratio transmittance with top level grid

// step how many voxels at once for indirect rays?
#define VOL_STEP_INDIRECT 1 // leave this at 1, it's broken otherwise.

typedef struct grid_cell_t
{
  float mu_max, mu_c, mu_r;
}
grid_cell_t;

typedef struct
{
  // scale * voxel = world space. scale is the length of one voxel edge
  float scale;
  // mean cosine
  float g;
  float emissive;
  // TODO: also need mu_t in spectral to convert density to color?
  float sigma_t;
  float sigma_s;

  // bounds for majorant and control variate:
  float max_density, avg_density;
  grid_cell_t *grid; // grid storing woodcock tracking related mu_*
  int grid_shift;    // ratio of downsampled grid size
  int grid_size[3];  // world space bounding box is the same as density, this is ss_aabb_density size
  int grid_coarse_size[3];  // super voxel resolution

  // bounding boxes in world space
  float ws_aabb_density[6];
  float ws_aabb_temperature[6];
  float ws_aabb_velocity[6];
  // bounding boxes in sample space
  int ss_aabb_density[6];
  int ss_aabb_temperature[6];
  int ss_aabb_velocity[6];
  int motion_blur;
  // total attenuation coefficient:
  openvdb::FloatGrid::Ptr density;
  // temperature/emission:
  openvdb::FloatGrid::Ptr temperature;
  // motion:
  openvdb::Vec3SGrid::Ptr velocity;
  int mshader;
}
medium_t;

#include "ovdblookup.h"
#include "woodcock.h"
#include "raymarch.h"

// switch on volume callbacks, to distinguish between hete and homo volumes:
extern "C" int volume_enabled(void *data)
{
  return 1;
}

namespace {

inline float blackbody_radiation(
    const float lambda,      // in [nm]
    const float temperature) // in [K]
{
  const double h = 6.6261e-34; // Planck's constant [J s]
  const double c = 2.9979e8; // speed of light [m/s]
  const double k = 1.3807e-23; // Boltzmann's constant [J/K]
  const double lambda_m = lambda*1e-9; // lambda [m]
  const double a = h*c/(k*temperature*lambda_m);
  const double lambda2 = lambda_m*lambda_m;
  const double lambda5 = lambda2*lambda_m*lambda2;
  const double b =(8.0*M_PI*h*c)/lambda5;
  return b/(exp(a)-1.0);
}

} // anonymous namespace


extern "C" float volume_transmittance(const path_t *p, int e, void *data)
{
#if VOL_METHOD == 1 // ray marching
  // pass in maximum length for clipping
  return volume_march<MODE_T>((medium_t *)data, p, e, p->e[e].dist);
#elif VOL_METHOD == 2 // woodcock tracking
  medium_t *s = (medium_t *)data;
  float tmin = 0.f, tmax = p->e[e].dist;
  // XXX might need to go to max velocity aabb.
  volume_get_interval_bounds(p->v[e-1].hit.x, p->e[e].omega, s->ws_aabb_density, &tmin, &tmax);
  if(tmin > tmax) return 1.0f;
  return woodcock_transmittance(
      s, p->v[e-1].hit.x, p->e[e].omega,
      s->sigma_t * s->avg_density, s->sigma_t * (s->max_density - s->avg_density), 0.0f, tmax, p->time, e == 1);
#elif VOL_METHOD == 3 // grid woodcock
  return woodcock_transmittance_grid((medium_t *)data, p->v[e-1].hit.x, p->e[e].omega, FLT_MAX, p->time, e==1);
#endif
}

extern "C" float volume_sample(path_t *p, int e, void *data)
{
  medium_t *s = (medium_t *)(data);
#if VOL_METHOD == 1 // ray marching
  return volume_march<MODE_SAMPLE>(s, p, e, pointsampler(p, s_dim_free_path));
#elif VOL_METHOD == 2 // global woodcock tracking
  float tmin = 0.f, tmax = FLT_MAX;
  // XXX might need to go to max velocity aabb.
  volume_get_interval_bounds(p->v[e-1].hit.x, p->e[e].omega, s->ws_aabb_density, &tmin, &tmax);
  if(tmin > tmax) return FLT_MAX;
  return woodcock_track(s, p->v[e-1].hit.x, p->e[e].omega,
      s->sigma_t * s->max_density, 0.0f, tmax, p->time, e==1);
#elif VOL_METHOD == 3 // grid woodcock
  return woodcock_track_grid(s, p->v[e-1].hit.x, p->e[e].omega, FLT_MAX, p->time, e==1);
#endif
}

extern "C" float volume_pdf_adjoint(const path_t *p, int e, void *data)
{
  return 1.0f;// XXX
  // quite the same as volume_pdf, except that we multiply mu_t from the other vertex.
  const float pdf = volume_march<MODE_PDF>((medium_t *)data, p, e, p->e[e].dist);
  if(primid_invalid(p->v[e-1].hit.prim) && !(p->v[e-1].flags & s_environment))
  {
    if(p->v[e-1].interior.mu_t < VOL_MU_T_CUTOFF) return 0.0f;
    return pdf * p->v[e-1].interior.mu_t;
  }
  return pdf;
}

// FIXME: pathspace actually calls volume_sample, volume_pdf, and volume_throughput in one _extend() call!!
// FIXME: unfortunately this is necessary, as we might have been clipped by geometry. might pay off to detect though.
// FIXME: actually i think it's just plain wrong in case geo is intersected before hand, need to check primid_valid in here!
extern "C" float volume_pdf(const path_t *p, int e, void *data)
{
  return 1.0f;// XXX
  const float pdf = volume_march<MODE_PDF>((medium_t *)data, p, e, p->e[e].dist);
  if(primid_invalid(p->v[e].hit.prim) && !(p->v[e].flags & s_environment))
  {
    if(p->v[e].interior.mu_t < VOL_MU_T_CUTOFF) return 0.0f;
    return pdf * p->v[e].interior.mu_t;
  }
  return pdf;
}

// called after init, in case we are somehow associated with a geo object.
// we interpret this to rescale to the object's bounding box.
extern "C" int uvn_init(const char *fname, uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = (medium_t *)self->data;
  s->mshader = self - rt.shader->shader;
  float aabb[6];
  prims_get_shape_aabb(rt.prims, shapeid, aabb);
  fprintf(stderr, "[medium_ovdb] captured shape[%d] (%s) bounding box %gx%gx%g\n",
      shapeid, fname,
      aabb[3]-aabb[0], aabb[4]-aabb[1], aabb[5]-aabb[2]);

  // need to copy transform, it might be shared and we don't want to apply our correction here twice.
  openvdb::Vec3d scale(
      (aabb[3]-aabb[0])/(s->ws_aabb_density[3]-s->ws_aabb_density[0]),
      (aabb[4]-aabb[1])/(s->ws_aabb_density[4]-s->ws_aabb_density[1]),
      (aabb[5]-aabb[2])/(s->ws_aabb_density[5]-s->ws_aabb_density[2]));
  openvdb::Vec3d trans0( -s->ws_aabb_density[0], -s->ws_aabb_density[1], -s->ws_aabb_density[2]);
  openvdb::Vec3d trans1( aabb[0], aabb[1], aabb[2]);

  openvdb::math::Transform *tr = new openvdb::math::Transform(s->density->constTransform());
  tr->postTranslate(trans0);
  tr->postScale(scale);
  tr->postTranslate(trans1);
  s->density->setTransform(openvdb::math::Transform::Ptr(tr));

  // recreate bounding boxes,
  // assume simple scale/translate, no rotations:
  openvdb::CoordBBox bbox = s->density->evalActiveVoxelBoundingBox();
  openvdb::Vec3d wm = s->density->constTransform().indexToWorld(bbox.min());
  openvdb::Vec3d wM = s->density->constTransform().indexToWorld(bbox.max());
  s->ws_aabb_density[0] = wm[0];
  s->ws_aabb_density[1] = wm[1];
  s->ws_aabb_density[2] = wm[2];
  s->ws_aabb_density[3] = wM[0];
  s->ws_aabb_density[4] = wM[1];
  s->ws_aabb_density[5] = wM[2];
  if(!s->density->transform().hasUniformScale())
    fprintf(stderr, "[medium_ovdb] rescaling object non-uniformly! this is probably not what you intended. (%g %g %g)\n", scale[0], scale[1], scale[2]);
  s->scale = s->density->transform().voxelSize()[0];

  if(s->emissive > 0.0f)
  {
    tr = new openvdb::math::Transform(s->temperature->constTransform());
    tr->postTranslate(trans0);
    tr->postScale(scale);
    tr->postTranslate(trans1);
    s->temperature->setTransform(openvdb::math::Transform::Ptr(tr));

    bbox = s->temperature->evalActiveVoxelBoundingBox();
    wm = s->temperature->constTransform().indexToWorld(bbox.min());
    wM = s->temperature->constTransform().indexToWorld(bbox.max());
    s->ws_aabb_temperature[0] = wm[0];
    s->ws_aabb_temperature[1] = wm[1];
    s->ws_aabb_temperature[2] = wm[2];
    s->ws_aabb_temperature[3] = wM[0];
    s->ws_aabb_temperature[4] = wM[1];
    s->ws_aabb_temperature[5] = wM[2];
  }

  if(s->motion_blur)
  {
    tr = new openvdb::math::Transform(s->velocity->constTransform());
    tr->postTranslate(trans0);
    tr->postScale(scale);
    tr->postTranslate(trans1);
    s->velocity->setTransform(openvdb::math::Transform::Ptr(tr));

    bbox = s->velocity->evalActiveVoxelBoundingBox();
    wm = s->velocity->constTransform().indexToWorld(bbox.min());
    wM = s->velocity->constTransform().indexToWorld(bbox.max());
    s->ws_aabb_velocity[0] = wm[0];
    s->ws_aabb_velocity[1] = wm[1];
    s->ws_aabb_velocity[2] = wm[2];
    s->ws_aabb_velocity[3] = wM[0];
    s->ws_aabb_velocity[4] = wM[1];
    s->ws_aabb_velocity[5] = wM[2];
  }
  return 1; // we want to discard the shape (only used for proxy aabb)
}

extern "C" float prepare(path_t *p, int v, void *data)
{ 
  medium_t *s = (medium_t *)data;
  // set volume properties on next segment
  p->e[v].vol.mean_cos = s->g;
  p->v[v].interior.mean_cos = s->g;
  p->v[v].interior.shader = s->mshader;
  p->v[v].material_modes = static_cast<vertex_scattermode_t>(s_volume | s_glossy);
  if(!(p->v[v].flags & s_environment) && primid_invalid(p->v[v].hit.prim))
  {
    // query density and derive mu_t, mu_s from that.
    openvdb::FloatGrid::ConstAccessor da = s->density->getConstAccessor();
    openvdb::Vec3SGrid::ConstAccessor va = s->velocity->getConstAccessor();
    const openvdb::math::Transform &dat = s->density->constTransform();
    const openvdb::math::Transform &vat = s->velocity->constTransform();
    const float density = volume_lookup_ws(p->v[v].hit.x, p->time, da, va, dat, vat, 0, s->motion_blur);
    p->v[v].interior.mu_t = density * s->sigma_t; //spectrum_rgb_to_p(p->lambda, s->mu_t);
    p->v[v].interior.mu_s = density * s->sigma_s;

    if(s->emissive > 0.0)
    {
      openvdb::FloatGrid::ConstAccessor ta = s->temperature->getConstAccessor();
      const openvdb::math::Transform &tat = s->temperature->constTransform();
      const float temperature = volume_lookup_ws(p->v[v].hit.x, p->time, ta, va, tat, vat, 0, s->motion_blur);
      if(temperature > 0.0)
      {
        p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_volume | s_emit | s_diffuse);
        p->v[v].shading.em = s->emissive * blackbody_radiation(p->lambda, temperature);
      }
    }
  }
  return 1.0f;
}

extern "C" float sample(path_t *p, void *data)
{
  const int v = p->length-1;
  if(p->v[v].interior.mu_t <= 0.0f) return 0.0f; // kill paths in the case that we sampled a position in vacuum
  // c++ sucks horse cock:
#ifdef __cplusplus
  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy | s_volume);
#else
  p->v[v].mode |= s_glossy | s_volume;
#endif
  const hit_t *hit = &p->v[v].hit;
  float out[3];
  sample_hg(p->v[v].interior.mean_cos, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), out, &p->v[v+1].pdf);
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = hit->n[k] * out[0] + hit->a[k] * out[1] + hit->b[k] * out[2];
  return p->v[v].interior.mu_s;
}

extern "C" float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_volume)) return 0.0f;
  return sample_eval_hg(p->v[v].interior.mean_cos, p->e[e1].omega, p->e[e2].omega);
}

extern "C" float brdf(path_t *p, int v, void *data)
{
  // c++, c++, ..
#ifdef __cplusplus
  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy | s_volume);
#else
  p->v[v].mode |= s_glossy | s_volume;
#endif
  if(p->v[v].interior.mu_t <= 0.0f) return 0.0f; // kill paths in the case that we sampled a position in vacuum
  return p->v[v].interior.mu_s * sample_eval_hg(p->v[v].interior.mean_cos, p->e[v].omega, p->e[v+1].omega);
}

extern "C" int init(FILE* f, void** data)
{
  medium_t *s = (medium_t *)malloc(sizeof(medium_t));
  memset(s, 0, sizeof(medium_t));
  *data = (void *)s;
  char filename[1024];
  s->sigma_t = 10.0f;
  s->sigma_s = 10.0f;
  s->emissive = 0.0;
  int i = fscanf(f, "%f %f %f %f %s", &s->g, &s->sigma_s, &s->sigma_t, &s->emissive, filename);
  if(i != 5)
  {
    fprintf(stderr, "[medium_ovdb] could not parse all arguments! expecting medium_ovdb <mean_cosine> <sigma_s> <sigma_t> <emissivity> <openvdb_filename>\n");
    s->g = 0.0f;
    return 1;
  }
  int res = fscanf(f, "%*[^\n]\n");

  openvdb::initialize(); // multiple calls are okay
  openvdb::io::File file(filename);
  try
  {
    file.open(); // read header only
  }
  catch (openvdb::IoError &e)
  {
    char fname[1024];
    snprintf(fname, 1024, "%s/%s", rt.searchpath, filename);
    file = openvdb::io::File(fname);
    file.open();
  }
  s->scale = 1.0;
  s->velocity = 0;
  for(openvdb::io::File::NameIterator it = file.beginName(); it != file.endName(); ++it)
  {
    if(it.gridName() == "density")
    {
      s->density = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(it.gridName()));
      openvdb::CoordBBox bbox = s->density->evalActiveVoxelBoundingBox();
      s->ss_aabb_density[0] = bbox.min()[0];
      s->ss_aabb_density[1] = bbox.min()[1];
      s->ss_aabb_density[2] = bbox.min()[2];
      s->ss_aabb_density[3] = bbox.max()[0];
      s->ss_aabb_density[4] = bbox.max()[1];
      s->ss_aabb_density[5] = bbox.max()[2];
      // assume simple scale/translate, no rotations:
      openvdb::Vec3d wm = s->density->constTransform().indexToWorld(bbox.min());
      openvdb::Vec3d wM = s->density->constTransform().indexToWorld(bbox.max());
      s->ws_aabb_density[0] = wm[0];
      s->ws_aabb_density[1] = wm[1];
      s->ws_aabb_density[2] = wm[2];
      s->ws_aabb_density[3] = wM[0];
      s->ws_aabb_density[4] = wM[1];
      s->ws_aabb_density[5] = wM[2];
      // s->density->transform().print(std::cout);
      assert(s->density->transform().hasUniformScale());
      s->scale = s->density->transform().voxelSize()[0];
      // fprintf(stderr, "scale = %g\n", s->scale);
    }
    if(it.gridName() == "temperature")
    {
      s->temperature = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(it.gridName()));
      openvdb::CoordBBox bbox = s->temperature->evalActiveVoxelBoundingBox();
      s->ss_aabb_temperature[0] = bbox.min()[0];
      s->ss_aabb_temperature[1] = bbox.min()[1];
      s->ss_aabb_temperature[2] = bbox.min()[2];
      s->ss_aabb_temperature[3] = bbox.max()[0];
      s->ss_aabb_temperature[4] = bbox.max()[1];
      s->ss_aabb_temperature[5] = bbox.max()[2];
      // assume simple scale/translate, no rotations:
      openvdb::Vec3d wm = s->temperature->constTransform().indexToWorld(bbox.min());
      openvdb::Vec3d wM = s->temperature->constTransform().indexToWorld(bbox.max());
      s->ws_aabb_temperature[0] = wm[0];
      s->ws_aabb_temperature[1] = wm[1];
      s->ws_aabb_temperature[2] = wm[2];
      s->ws_aabb_temperature[3] = wM[0];
      s->ws_aabb_temperature[4] = wM[1];
      s->ws_aabb_temperature[5] = wM[2];
      // s->temperature->transform().print(std::cout);
    }
    if(it.gridName() == "v" || it.gridName() == "velocity")
    {
      s->velocity = openvdb::gridPtrCast<openvdb::Vec3SGrid>(file.readGrid(it.gridName()));
      openvdb::CoordBBox bbox = s->velocity->evalActiveVoxelBoundingBox();
      s->ss_aabb_velocity[0] = bbox.min()[0];
      s->ss_aabb_velocity[1] = bbox.min()[1];
      s->ss_aabb_velocity[2] = bbox.min()[2];
      s->ss_aabb_velocity[3] = bbox.max()[0];
      s->ss_aabb_velocity[4] = bbox.max()[1];
      s->ss_aabb_velocity[5] = bbox.max()[2];
      // assume simple scale/translate, no rotations:
      openvdb::Vec3d wm = s->velocity->constTransform().indexToWorld(bbox.min());
      openvdb::Vec3d wM = s->velocity->constTransform().indexToWorld(bbox.max());
      s->ws_aabb_velocity[0] = wm[0];
      s->ws_aabb_velocity[1] = wm[1];
      s->ws_aabb_velocity[2] = wm[2];
      s->ws_aabb_velocity[3] = wM[0];
      s->ws_aabb_velocity[4] = wM[1];
      s->ws_aabb_velocity[5] = wM[2];
      // s->velocity->transform().print(std::cout);
    }
    // TODO: what is `fuel' ? float grid with slightly smaller grid dimensions?
  }
  file.close();
  if(!s->velocity)
  {
    s->velocity = openvdb::Vec3SGrid::Ptr(new openvdb::Vec3SGrid());
    s->motion_blur = 0;
  }
  else s->motion_blur = 1;
  // XXX
  s->motion_blur = 0;

  // now get some composite extinction coefficients:
  double start = common_time_wallclock();
  // downsampling ratio. choose this big in almost homogeneous volumes.
  // will result in low performance (samples/second) if too small (voxels too small),
  // and will result in poor variance if too large (contents too heterogeneous).
  s->grid_shift = 7;
  const int wd = (1<<s->grid_shift);
  const int size[3] = {
    s->ss_aabb_density[3]-s->ss_aabb_density[0],
    s->ss_aabb_density[4]-s->ss_aabb_density[1],
    s->ss_aabb_density[5]-s->ss_aabb_density[2]};
  for(int k=0;k<3;k++) s->grid_size[k] = size[k];
  for(int k=0;k<3;k++) s->grid_coarse_size[k] = (size[k]-1)/wd + 1;
  fprintf(stderr, "[medium ovdb] computing control variate for extinction on a %dx%dx%d - %dx%dx%d grid\n", s->grid_coarse_size[0], s->grid_coarse_size[1], s->grid_coarse_size[2], size[0], size[1], size[2]);
  s->grid = (grid_cell_t *)common_alloc(16, s->grid_coarse_size[0]*s->grid_coarse_size[1]*s->grid_coarse_size[2]*sizeof(grid_cell_t));
  memset(s->grid, 0, s->grid_coarse_size[0]*s->grid_coarse_size[1]*s->grid_coarse_size[2]*sizeof(grid_cell_t));

  openvdb::FloatGrid::ConstAccessor a = s->density->getConstAccessor();
  double max_density = 0.0f, avg_density = 0.0f;
  size_t zeroes = 0;
  const float vs[3] = {
    (s->ws_aabb_density[3]-s->ws_aabb_density[0])/s->grid_size[0],
    (s->ws_aabb_density[4]-s->ws_aabb_density[1])/s->grid_size[1],
    (s->ws_aabb_density[5]-s->ws_aabb_density[2])/s->grid_size[2]};
  const float D = sqrtf(dotproduct(vs,vs))*(float)wd; // voxel diagonal

  for(int k=0;k<s->grid_coarse_size[2];k++)
  for(int j=0;j<s->grid_coarse_size[1];j++)
  for(int i=0;i<s->grid_coarse_size[0];i++)
  {
    float mu_max = 0.0f, mu_avg = 0.0f, mu_min = FLT_MAX;
    int iii, jjj, kkk, cnt = 0;
    for(int kk=0;kk<wd && (kkk = k*wd+kk+s->ss_aabb_density[2]) < s->ss_aabb_density[5];kk++)
    for(int jj=0;jj<wd && (jjj = j*wd+jj+s->ss_aabb_density[1]) < s->ss_aabb_density[4];jj++)
    for(int ii=0;ii<wd && (iii = i*wd+ii+s->ss_aabb_density[0]) < s->ss_aabb_density[3];ii++)
    {
      cnt ++;
      float density = a.getValue(openvdb::Coord(iii,jjj,kkk));
      if(density < VOL_MU_T_CUTOFF) density = 0.0f;
      max_density = fmaxf(density, max_density);
      avg_density += density;
      mu_max = fmaxf(density, mu_max);
      mu_min = fminf(density, mu_min);
      mu_avg += density;
    }
    mu_avg /= (float)cnt;
    if(mu_avg == 0.0) zeroes++;

    // apply cross sections to density to obtain real extinctions:
    mu_avg *= s->sigma_t;
    mu_min *= s->sigma_t;
    mu_max *= s->sigma_t;
    // apply eq(6) from [Novak et al. 2014]. does not fix my problems.
    const float gamma = 2.0f;
    float mu_c = mu_avg;
    float mu_r = fmaxf(mu_c - mu_min, mu_max - mu_c);
    mu_c = mu_min + mu_r * (powf(gamma, 1.0f/(D*mu_r)) - 1.0f);
    mu_c = fminf(fmaxf(mu_min, mu_c), mu_avg);
    // fprintf(stderr, "vox %d %d %d with [%f %f %f] %f\n", i, j, k, mu_min, mu_avg, mu_max, mu_c);

    // store mu_r and mu_c and mu_max
    s->grid[i + s->grid_coarse_size[0]*(j + k*s->grid_coarse_size[1])].mu_r = fmaxf(mu_c - mu_min, mu_max - mu_c);
    s->grid[i + s->grid_coarse_size[0]*(j + k*s->grid_coarse_size[1])].mu_c = mu_c;
    s->grid[i + s->grid_coarse_size[0]*(j + k*s->grid_coarse_size[1])].mu_max = mu_max;
  }
  avg_density /= (double)size[0]*(double)size[1]*(double)size[2];
  s->max_density = max_density;
  s->avg_density = avg_density;
  double end = common_time_wallclock();
  fprintf(stderr, "[medium ovdb] %zu empty, average (%g) and max density (%g). took %.03fs.\n", zeroes, avg_density, max_density, end-start);

  return res == -1;
}

extern "C" void cleanup(void *data)
{
  // ovdb is all in shared pointers, so we're good.
  medium_t *s = reinterpret_cast<medium_t *>(data);
  free(s->grid);
}

#undef VOL_METHOD
#undef VOL_MU_T_CUTOFF
#undef MODE_T
#undef MODE_SAMPLE
#undef MODE_PDF

