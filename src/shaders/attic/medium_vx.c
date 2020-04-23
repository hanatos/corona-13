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
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

// heterogeneous medium, using a stupid full voxel grid as backend.
// 
// this doesn't perform woodcock tracking but works on the assumption that
// voxels represent piecewise homogeneous media. conceptually, this code
// constructs a cdf for every voxel along a ray, samples from that, and
// then samples a distance within this homogeneous block. in practice, the
// random number is chosen up front, so the cdf is never computed as whole.

#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

// volume positions with lower mu_t concentration are considered empty:
#define VOL_MU_T_CUTOFF 1e-5f


typedef struct vx_header_t
{
  int v, x;            // file magic number
  int32_t size[3];     // voxel counts for each dimension
  float aabb[6];       // world space bounding box
  float *data;         // mmapped data pointer
}
vx_header_t;

typedef struct
{
  // scale * voxel = world space. scale is the length of one voxel edge
  float scale;
  // mean cosine
  float g;
  int emissive;
  // TODO: also need mu_t in spectral to convert density to color?
  float sigma_t;
  float sigma_s;
  float brightness;
  int mshader;

  int vx_file;
  size_t vx_data_size;
  void *vx_data;
  vx_header_t vx;
}
medium_t;

static inline float blackbody_radiation(
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

// returns the start/end interval in (fractional) voxel index space coordinates,
// given the path edge p->e[e] and the world space voxel bounding box aabb
static inline float volume_get_interval(
    const vx_header_t *vx,
    const path_t *p,
    const int e,
    float tmax,
    float *start,
    float *end)
{
  // intersect aabb
  float tmin = 0.f;
  for(int k=0;k<3;k++)
  {
    const float t0 = (vx->aabb[k+3] - p->v[e-1].hit.x[k])/p->e[e].omega[k];
    const float t1 = (vx->aabb[k]   - p->v[e-1].hit.x[k])/p->e[e].omega[k];
    if(t0 <= t1)
    {
      tmin = t0 > tmin ? t0 : tmin;
      tmax = t1 < tmax ? t1 : tmax;
    }
    else
    {
      tmin = t1 > tmin ? t1 : tmin;
      tmax = t0 < tmax ? t0 : tmax;
    }
  }
  // can only work if start != end because only then we get a valid step direction:
  if(tmin < tmax)
  {
    // compute world space start/end point. already subtract aabb min:
    for(int k=0;k<3;k++) start[k] = p->v[e-1].hit.x[k] + p->e[e].omega[k] * tmin;
    for(int k=0;k<3;k++) end  [k] = p->v[e-1].hit.x[k] + p->e[e].omega[k] * tmax;
    // convert to sample space (in voxel grid):
    for(int k=0;k<3;k++) start[k] = (start[k] - vx->aabb[k])*vx->size[k]/(vx->aabb[3+k] - vx->aabb[k]);
    for(int k=0;k<3;k++) end  [k] = (end  [k] - vx->aabb[k])*vx->size[k]/(vx->aabb[3+k] - vx->aabb[k]);
    assert(tmin >= 0);
    return tmin;
  }
  return -1;
}

static inline size_t volume_index(const vx_header_t *vx, const float *pos)
{
  ssize_t ind = 5*(vx->size[0]*(vx->size[1]*(int64_t)pos[2] + (int64_t)pos[1]) + (int64_t)pos[0]);
  assert(ind >= 0);
  assert(ind < (int64_t)vx->size[2]*vx->size[1]*vx->size[0]*5);
  return ind;
}

#define DENSITY 3
#define TEMPERATURE 4
static inline float volume_lookup_ss(
    const vx_header_t *vx,
    const int which,                           // 3 is density, 4 is temperature. 0-2 are velocity.
    const float *pos,                          // lookup coordinates in local index space of accessor a
    const float time,                          // path's time for motion blur
    const int eyeray)                          // flag to signify high quality interpolation or not
{
  // unfiltered:
  for(int k=0;k<3;k++) if(pos[k] < 0 || pos[k] >= vx->size[k]) return 0.0f;
  const float *v = vx->data + volume_index(vx, pos);

  float advected[3];
  for(int k=0;k<3;k++) advected[k] = pos[k] + v[k]*(time - .5);
  for(int k=0;k<3;k++) if(advected[k] < 0 || advected[k] >= vx->size[k]) return 0.0f;
  return vx->data[volume_index(vx, advected) + which];
}

static inline float volume_lookup_ws(
    const vx_header_t *vx,
    const int which,
    const float *pos,                          // lookup coordinates in world space
    const float time,                          // path's time for motion blur
    const int eyeray)                          // flag to signify high quality interpolation or not
{
  float sspos[3];
  for(int k=0;k<3;k++) sspos[k] = (pos[k] - vx->aabb[k])*vx->size[k]/(vx->aabb[3+k] - vx->aabb[k]);
  return volume_lookup_ss(vx, which, sspos, time, eyeray);
}

#define MODE_T 0
#define MODE_SAMPLE 1
#define MODE_PDF 2
// MODE:
// T      - compute and return transmittance, use rdist as tmax
// SAMPLE - use rdist as random to compute scattering distance 
// PDF    - use rdist as distance to compute pdf
static inline float volume_march(
    const medium_t *s,    // medium shader struct
    const path_t *p,      // light transport path
    const int e,          // edge number to march
    const int MODE,       // transmission, sample, or pdf?
    const float rdist)    // random number to sample the bin or distance
{
  float transmittance = 0.0f;
  float cdf = 0.0f;
  int i = 0;

  float start[3], end[3];
  const float tmin = volume_get_interval(
      &s->vx, p, e,
      MODE==MODE_T?rdist:FLT_MAX,
      start, end);
  if(MODE==MODE_SAMPLE && tmin < 0.0f) return FLT_MAX;
  if(MODE==MODE_PDF    && tmin < 0.0f) return 1.0;
  if(MODE==MODE_T      && tmin < 0.0f) return 1.0;

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
  const float dist = s->scale * sqrtf(dotproduct(step,step));
  const int pdf_bin = (MODE==MODE_PDF) ? (rdist == FLT_MAX ? 1<<30 : (rdist-tmin)/dist): 0.0;
  assert(MODE!=MODE_PDF || pdf_bin >= 0);
  float curr[3] = {start[0], start[1], start[2]};
  int eyeray = (p->v[e-1].mode & s_sensor); // fails for light tracing :(
  while(1)
  {
    // quit loop:
    for(int k=0;k<3;k++)
    {
      if(step[k] > 0 && (int)curr[k] >= s->vx.size[k]) goto break_outer;
      if(step[k] < 0 && (int)curr[k] < 0) goto break_outer;
    }
    if(step[dim] > 0 && (int)curr[dim] > (int)end[dim]) goto break_outer;
    if(step[dim] < 0 && (int)curr[dim] < (int)end[dim]) goto break_outer;

    // switch off filtering after 50% of the light transmittance is gone:
    if(eyeray && (transmittance > -logf(.5f)/dist)) eyeray = 0;
    const float mu_t = s->sigma_t * volume_lookup_ss(&s->vx, DENSITY, curr, p->time, eyeray);

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
    
    for(int k=0;k<3;k++) curr[k] += step[k];
  }
break_outer:
  // exhausted voxels while searching for cdf > rand
  if(MODE == MODE_SAMPLE)
    return FLT_MAX;
  // pdf to go through the whole volume is transmittance, too (MODE_T and MODE_PDF)
  return expf(-transmittance*dist);
}



// switch on volume callbacks, to distinguish between hete and homo volumes:
int volume_enabled(void *data)
{
  return 1;
}

float volume_transmittance(const path_t *p, int e, void *data)
{
  // pass in maximum length for clipping
  return volume_march((medium_t *)data, p, e, MODE_T, p->e[e].dist);
}

float volume_sample(path_t *p, int e, void *data)
{
  return volume_march((medium_t *)data, p, e, MODE_SAMPLE, pointsampler(p, s_dim_free_path));
}

float volume_pdf_adjoint(const path_t *p, int e, void *data)
{
  // quite the same as volume_pdf, except that we multiply mu_t from the other vertex.
  const float pdf = volume_march((medium_t *)data, p, e, MODE_PDF, p->e[e].dist);
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
float volume_pdf(const path_t *p, int e, void *data)
{
  const float pdf = volume_march((medium_t *)data, p, e, MODE_PDF, p->e[e].dist);
  if(primid_invalid(p->v[e].hit.prim) && !(p->v[e].flags & s_environment))
  {
    if(p->v[e].interior.mu_t < VOL_MU_T_CUTOFF) return 0.0f;
    return pdf * p->v[e].interior.mu_t;
  }
  return pdf;
}

float volume_throughput(const path_t *p, int e, void *data)
{
  return 1.0f;
}


int shape_init(uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = (medium_t *)self->data;
  s->mshader = self - rt->shader->shader;
  return 0;
}

float prepare(path_t *p, int v, void *data)
{ 
  medium_t *s = (medium_t *)data;
  // set volume properties on next segment
  p->e[v].vol.mean_cos = s->g;
  p->v[v].interior.mean_cos = s->g;
  p->v[v].interior.shader = s->mshader;
  p->v[v].material_modes = s_volume | s_glossy;
  if(!(p->v[v].flags & s_environment) && primid_invalid(p->v[v].hit.prim))
  {
    // query density and derive mu_t, mu_s from that.
    const float density = volume_lookup_ws(&s->vx, DENSITY, p->v[v].hit.x, p->time, 0);
    const float temperature = volume_lookup_ws(&s->vx, TEMPERATURE, p->v[v].hit.x, p->time, 0);
    // XXX mu_t and rv feel like they should be on the vertex!
    p->e[v].vol.mu_t =
    p->v[v].interior.mu_t = density * s->sigma_t; //spectrum_rgb_to_p(p->lambda, s->mu_t);
    // XXX sort those terms better! 
    p->v[v].interior.rv =
    p->e[v].vol.rv = s->sigma_s/s->sigma_t; // density cancels out as we're doing perfect importance sampling.
    if(s->emissive && temperature > 0.0)
    {
      p->v[v].mode |= s_volume | s_emit | s_diffuse;
      p->v[v].shading.em = s->brightness * blackbody_radiation(p->lambda, temperature);
    }
  }
  return 1.0f;
}

float sample(path_t *p, void *data)
{
  const int v = p->length-1;
  if(p->v[v].interior.mu_t < VOL_MU_T_CUTOFF) return 0.0f; // kill paths in the case that we sampled a position in vacuum
  p->v[v].mode |= s_glossy | s_volume;
  const hit_t *hit = &p->v[v].hit;
  float out[3];
  sample_hg(p->v[v].interior.mean_cos, pointsampler(p, s_dim_omega_x), pointsampler(p, s_dim_omega_y), out, &p->v[v+1].pdf);
  for(int k=0;k<3;k++) p->e[v+1].omega[k] = hit->n[k] * out[0] + hit->a[k] * out[1] + hit->b[k] * out[2];
  return p->v[v].interior.rv;
}

float pdf(path_t *p, int e1, int v, int e2, void *data)
{
  if(!(p->v[v].mode & s_volume)) return 0.0f;
  return sample_eval_hg(p->v[v].interior.mean_cos, p->e[e1].omega, p->e[e2].omega);
}

float brdf(path_t *p, int v, void *data)
{
  p->v[v].mode |= s_glossy | s_volume;
  if(p->v[v].interior.mu_t < VOL_MU_T_CUTOFF) return 0.0f; // kill paths in the case that we sampled a position in vacuum
  return p->v[v].interior.rv * sample_eval_hg(p->v[v].interior.mean_cos, p->e[v].omega, p->e[v+1].omega);
}

int init(FILE* f, void** data)
{
  medium_t *s = (medium_t *)malloc(sizeof(medium_t));
  memset(s, 0, sizeof(medium_t));
  *data = (void *)s;
  char filename[1024];
  s->sigma_t = 10.0f;
  s->sigma_s = 10.0f;
  s->emissive = 0;
  int i = fscanf(f, "%f %f %f %d %s", &s->g, &s->sigma_s, &s->sigma_t, &s->emissive, filename);
  if(i != 5)
  {
    fprintf(stderr, "[medium_vx] could not parse all arguments! expecting medium_vx <mean_cosine> <sigma_s> <sigma_t> <emit?> <vx_filename>\n");
    s->g = 0.0f;
    return 1;
  }
  int res = fscanf(f, "%*[^\n]\n");

  s->brightness = 1.f;

  int fd = open(filename, O_RDONLY);
  if(fd < 0)
  {
    char fname[1024];
    snprintf(fname, 1024, "%s/%s", rt.searchpath, filename);
    fd = open(fname, O_RDONLY);
  }
  if(fd < 0)
  {
    fprintf(stderr, "[medium_vx] could not open `%s'!\n", filename);
    return 1;
  }

  size_t data_size = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, SEEK_SET);
  s->vx_data = mmap(0, data_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
  memcpy(&s->vx, s->vx_data, sizeof(vx_header_t));
  s->vx.data = ((float *)s->vx_data) + 11;
  s->vx_file = fd;
  s->vx_data_size = data_size;

  // needs to be isotropic, or else we die.
  s->scale = (s->vx.aabb[3]-s->vx.aabb[0])/s->vx.size[0];

  return res == -1;
}

void cleanup(void *data)
{
  medium_t *s = (medium_t *)data;
  if(s->vx_data) munmap(s->vx_data, s->vx_data_size);
  if(s->vx_file > 2) close(s->vx_file);
  s->vx_file = -1;
  free(s);
}

#undef VOL_MU_T_CUTOFF
#undef DENSITY
#undef TEMPERATURE
#undef MODE_T
#undef MODE_SAMPLE
#undef MODE_PDF

