#pragma once

#include "pathspace.h"
#include "rgb2spec.h"
#include "shader.h"

// texture slot, i.e. where to put it on the path struct
typedef enum tex_slot_t
{
  s_slot_diffuse,
  s_slot_specular,
  s_slot_emission,
  s_slot_volume,
  s_slot_glossy,
  s_slot_roughness,
  s_slot_transmit_to_eye,
  s_slot_diffgeo_n,
  s_slot_diffgeo_dndu,
  s_slot_diffgeo_dndv,
  s_slot_unused,          // set to this if you don't want your data to be filled
  // ?? how to store dpdu/dpdv as per-vertex data? orthogonalise dpdu/dpdv around new n?
  s_slot_count
}
tex_slot_t;

static inline tex_slot_t tex_parse_slot(const char c)
{
  // TODO: multi-res diffgeo!
  // TODO: normals! dpdu/dpdv gn n
  if(c == 'd') return s_slot_diffuse;
  if(c == 's') return s_slot_specular;
  if(c == 'e') return s_slot_emission;
  if(c == 'v') return s_slot_volume;
  if(c == 'g') return s_slot_glossy;
  if(c == 'r') return s_slot_roughness;
  if(c == 't') return s_slot_transmit_to_eye;
  if(c == 'x') return s_slot_unused;
  return s_slot_count;
}

static inline void tex_set_slot(path_t *p, int v, tex_slot_t s, mf_t data)
{
  switch(s)
  {
    case s_slot_diffuse:   p->v[v].shading.rd        = data; return;
    case s_slot_specular:  p->v[v].shading.rs        = data; return;
    case s_slot_glossy:    p->v[v].shading.rg        = data; return;
    case s_slot_volume:
                           // for homogeneous interior media, we set the lobe modes here, the medium
                           // shader's prepare() is only called at the surface!
                           p->v[v].material_modes    = s_volume | s_glossy;
                           p->v[v].interior.mu_s     = data;
                           p->v[v].interior.mu_t     = mf_set1(1.0);
                           return;
    case s_slot_roughness: p->v[v].shading.roughness = mf(data, 0); return;
    case s_slot_emission:  p->v[v].shading.em        = data; return;
    case s_slot_transmit_to_eye:
      // set glossy colour only when transmitting towards the eye.
      if((p->v[0].mode & s_sensor) && !(p->v[v].flags & s_inside))
        p->v[v].shading.rg = data;
      else if((p->v[0].mode & s_emit) && (p->v[v].flags & s_inside))
        p->v[v].shading.rg = data;
      return;
    default:
      return;
  }
}

static inline void tex_set_slot_coeff(
    path_t *p, int v, tex_slot_t s, float mul, const float *coeff)
{
  mf_t val = mf_mul(mf_set1(mul), mf_rgb2spec(coeff, p->lambda));
  switch(s)
  {
    case s_slot_transmit_to_eye:
    case s_slot_diffuse:
    case s_slot_specular:
    case s_slot_glossy:
    case s_slot_volume:    tex_set_slot(p, v, s, mf_clamp(val, 0.0f, 1.0f)); return;
    case s_slot_roughness: tex_set_slot(p, v, s, mf_clamp(val, 0.0f, 1.0f)); return;
    case s_slot_emission:  tex_set_slot(p, v, s, val); return;
    default: return;
  }
}
