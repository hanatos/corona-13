// heterogeneous medium, wrapping around tuv container

extern "C"
{
#include "corona_common.h"
#include "shader.h"
#include "sampler_common.h"
#include "spectrum_common.h"
#include "prims.h"
}

#include "vol/Raytrace.hh"
#include "vol/VolumeStorage.hh"
#include "vol/VolumeTree.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

  vol::VolumeStorage *storage;
  vol::VolumeTree tree;
  int mshader;
}
medium_t;

// switch on volume callbacks, to distinguish between hete and homo volumes:
extern "C" int volume_enabled(void *data)
{
  return 1;
}

static inline float blackbody_radiation(
    const float lambda,      // in [nm]
    const float temperature) // in [K]
{
  const double h = 6.62606957e-34; // Planck's constant [J s]
  const double c = 299792458.0;    // speed of light [m/s]
  const double k = 1.3807e-23;     // Boltzmann's constant [J/K]
  const double lambda_m = lambda*1e-9; // lambda [m]
  const double lambda2 = lambda_m*lambda_m;
  const double lambda5 = lambda2*lambda_m*lambda2;
  const double c1 = 2. * h * c * c / lambda5;
  const double c2 = h * c / (lambda_m * temperature * k);
  // convert to spectral radiance in [W/m^2 / sr / nm]
  return c1 / (exp(c2)-1.0) * 1e-9;
}

extern "C" float volume_transmittance(const path_t *p, int e, void *data)
{
  medium_t *s = (medium_t *)data;
  return vol::EvalTransmittance(s->tree, p->v[e-1].hit.x, p->e[e].omega, p->time, s->sigma_t, 0.0f, p->e[e].dist);
}

extern "C" float volume_sample(path_t *p, int e, void *data)
{
  medium_t *s = (medium_t *)data;
  return vol::SampleFreePath(s->tree, p->v[e-1].hit.x, p->e[e].omega, p->time, s->sigma_t, 0.0f, FLT_MAX, pointsampler(p, s_dim_free_path));
}

extern "C" float volume_pdf_adjoint(const path_t *p, int e, void *data)
{
  // TODO: this is something like transmittance, or transmittance * mu_t on some side?
  return 1.0f;// XXX
}

// FIXME: pathspace actually calls volume_sample, volume_pdf, and volume_throughput in one _extend() call!!
// FIXME: unfortunately this is necessary, as we might have been clipped by geometry. might pay off to detect though.
// FIXME: actually i think it's just plain wrong in case geo is intersected before hand, need to check primid_valid in here!
extern "C" float volume_pdf(const path_t *p, int e, void *data)
{
  return 1.0f;// XXX
}

extern "C" float volume_throughput(const path_t *p, int e, void *data)
{
  if(p->v[e].material_modes & s_volume) return 1.0f/p->v[e].interior.mu_t;
  return 1.0f;
}


// called after init, in case we are somehow associated with a geo object.
// we interpret this to rescale to the object's bounding box.
extern "C" int uvn_init(const char *fname, uint32_t shapeid, shader_so_t *self)
{
  medium_t *s = (medium_t *)self->data;
  s->mshader = self - rt->shader->shader;
#if 0 // TODO: add support for this
  float aabb[6];
  prims_get_shape_aabb(rt->prims, shapeid, aabb);
  fprintf(stderr, "[medium_hete] captured shape[%d] (%s) bounding box %gx%gx%g\n",
      shapeid, fname,
      aabb[3]-aabb[0], aabb[4]-aabb[1], aabb[5]-aabb[2]);

  // TODO: allow rotations?
  // adjust voxel size to new scale, only use first dimension (rest isotropic scaling)
  const float scale = (aabb[3]-aabb[0])/(s->tree->content_box[3]-s->tree->content_box[0]);
  s->tree->voxel_size = s->tree->header->voxel_size * scale;
  float loc[3];
  for(int k=0;k<3;k++) loc[k] = aabb[k] - s->tree->aabb[k];
  for(int k=0;k<6;k++) s->tree->content_box[k] *= scale;
  for(int k=0;k<6;k++) s->tree->aabb[k] *= scale;
  vol_create_transform(s->tree, loc, 0);
  return 1; // we want to discard the shape (only used for proxy aabb)
#endif
  return 0;
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
    vol::Payload payload;
    vol::Lookup(s->tree, p->v[v].hit.x, p->time, &payload);

    // XXX mu_t and mu_s feel like they should be on the vertex!
    p->e[v].vol.mu_t =
    p->v[v].interior.mu_t = payload.m_density * s->sigma_t; //spectrum_rgb_to_p(p->lambda, s->mu_t);
    // XXX sort those terms better! 
    p->v[v].interior.mu_s =
    p->e[v].vol.mu_s = payload.m_density * s->sigma_s;

    if(s->emissive > 0.0)
    {
      if(payload.m_temperature > 0.0)
      {
        p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_volume | s_emit | s_diffuse);
        p->v[v].shading.em = s->emissive * blackbody_radiation(p->lambda, payload.m_temperature);
      }
    }
  }
  return 1.0f;
}

extern "C" float sample(path_t *p, void *data)
{
  const int v = p->length-1;
  if(p->v[v].interior.mu_t <= 0.0f) return 0.0f; // kill paths in the case that we sampled a position in vacuum
  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy | s_volume);
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
  p->v[v].mode = static_cast<vertex_scattermode_t>(p->v[v].mode | s_glossy | s_volume);
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
    fprintf(stderr, "[medium_tuv] could not parse all arguments! expecting medium_tuv <mean_cosine> <sigma_s> <sigma_t> <emissivity> <tree_filename>\n");
    s->g = 0.0f;
    return 1;
  }
  int res = fscanf(f, "%*[^\n]\n");

  s->storage = vol::VolumeStorage::Open('r', filename);
  if(!s->storage)
  {
    char fname[1024];
    snprintf(fname, 1024, "%s/%s", rt.searchpath, filename);
    s->storage = vol::VolumeStorage::Open('r', fname);
  }
  if(!s->storage)
  {
    fprintf(stderr, "[medium_tuv] failed to load tuv file `%s'\n", filename);
    return 1;
  }
  fprintf(stdout, "[medium_tuv] loaded tuvfile magic %lu version %lu\n", s->storage->m_header->magic, s->storage->m_header->version);

  // get kd tree:
  s->storage->Map(s->tree);
  s->scale = 1.0;

  return res == -1;
}

extern "C" void cleanup(void *data)
{
  medium_t *s = (medium_t *)data;

  // TODO: destroy storage and tree
  s->storage->Close();
}
