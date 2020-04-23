#include "lights.h"
#include "light_hierarchy.h"

typedef struct lights_t
{
  lh_t lh;
  float p_geo, p_sky, p_vol;
  shader_so_t *vol;
}
lights_t;

// TODO: init this in better ways
void lights_pdf_type(const path_t *p, const int v, mf_t p_egv[3])
{
  p_egv[0] = mf_set1(rt.lights->p_sky);
  p_egv[1] = mf_set1(rt.lights->p_geo);
  p_egv[2] = mf_set1(rt.lights->p_vol);
}

// global initialisation
lights_t *lights_init()
{
  lights_t *s = malloc(sizeof(*s));
  memset(s, 0, sizeof(*s));
  s->p_geo = s->p_sky = s->p_vol = 0.0f;
  lh_init(&s->lh);
  return s;
}

// ..and cleanup
void lights_cleanup(lights_t *s)
{
  lh_cleanup(&s->lh);
  free(s);
}

// initialise the give shapeid as light source, use L as average emitted
// radiance for importance sampling
void lights_init_light(const uint32_t shapeid, const float L)
{
  uint64_t off = rt.lights->lh.num_prims;
  rt.lights->lh.num_prims += rt.prims->shape[shapeid].num_prims;
  if(rt.lights->lh.max_num_prims < rt.lights->lh.num_prims)
  {
    rt.lights->lh.prims = realloc(rt.lights->lh.primid, sizeof(lh_prim_t)*rt.lights->lh.num_prims);
    rt.lights->lh.max_num_prims = rt.lights->lh.num_prims;
  }
  for(int k=0;k<rt.prims->shape[shapeid].num_prims;k++)
  {
    primid_t pi = rt.prims->shape[shapeid].primid[k];
    pi.shapeid = shapeid; // can't touch the mmapped primid on the shape, so we need to update here
    if(pi.vcnt == 3 || pi.vcnt == 4)
    { // discard non tri/quad light sources for now
      rt.lights->lh.prims[off].primid = pi;
      rt.lights->lh.prims[off].energy = L * prims_get_area(rt.prims, pi);
      float c[3] = {0.0f};
      float aabb[6] = {FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};
      for(int i=0;i<pi.vcnt;i++)
      {
        const float *v = geo_get_vertex(rt.prims, pi, 0);
        for(int j=0;j<3;j++)
        {
          aabb[j]   = MIN(aabb[j],   v[j]);
          aabb[3+j] = MAX(aabb[3+j], v[j]);
          c[j] += v[j];
        }
      }
      for(int i=0;i<pi.vcnt;i++) c[j] /= pi.vcnt;
      lh_cone_t cone = lh_triangle_to_cone(
          geo_get_vertex(rt.prims, pi, 0),
          geo_get_vertex(rt.prims, pi, 1),
          geo_get_vertex(rt.prims, pi, 2));
      if(pi.vcnt == 4)
      {
        lh_cone_t cone1 = lh_triangle_to_cone(
            geo_get_vertex(rt.prims, pi, 0),
            geo_get_vertex(rt.prims, pi, 2),
            geo_get_vertex(rt.prims, pi, 3));
        cone = lh_cone_union(cone, cone1);
      }
      rt.lights->lh.prims[off].cone = cone;
      memcpy(rt.lights->lh.prims[off].c, c, sizeof(float)*3);
      memcpy(rt.lights->lh.prims[off].aabb, aabb, sizeof(float)*6);
      for(int j=0;j<3;j++)
      { // adjust root lh bounding box
        rt.lights->lh.aabb[j]   = MIN(rt.lights->lh.aabb[j], aabb[j]);
        rt.lights->lh.aabb[3+j] = MAX(rt.lights->lh.aabb[3+j], aabb[3+j]);
      }
      off++;
    }
  }
}

// called once before every progression
void lights_prepare_frame()
{
  if(!rt.lights->lh.nodes)
  { // not initialised yet
    lh_build_binned(&rt.lights->lh);
  }
}


// computes the spectral pdf of next event estimation for p->v[v]
mf_t lights_pdf_next_event(const path_t *p, int v);

// returns throughput into direction of segment, no visibility is checked here.
// will construct vertex p->v[p->length] and not increment length
mf_t lights_sample_next_event(path_t *p);

// sample a light ray at p->v[0] and p->e[1]
mf_t lights_sample(path_t *p);

// returns the pdf of sampling the first or last segment
// of the path (depending on v) via lights_sample()
mf_t lights_pdf(const path_t *p, int v);

// evaluate emitted radiance of given path vertex v,
// in direction towards the sensor.
mf_t lights_eval_vertex(path_t *path, int v);




// volume sampling wrappers

// attach a volume emitter to the global list
void lights_init_volume_light(shader_so_t *vol)
{
  rt.lights->vol = vol;
}

mf_t light_volume_pdf_nee(const path_t *p, int v)
{
  return rt.lights->vol->volume_pdf_nee(p, v, rt.lights->vol->data);
}

mf_t light_volume_pdf_fnee_direction(const path_t *p, int v)
{
  return rt.lights->vol->volume_pdf_fnee_direction(p, v, rt.lights->vol->data);
}

mf_t light_volume_pdf(const path_t *p, int v)
{
  return rt.lights->vol->volume_pdf(p, v, rt.lights->vol->data);
}

mf_t light_volume_sample_nee(path_t *p, int v)
{
  return rt.lights->vol->volume_sample_nee(p, v, rt.lights->vol->data);
}

void light_volume_sample_fnee_direction(path_t *p, int v)
{
  rt.lights->vol->volume_sample_fnee_direction(p, v, rt.lights->vol->data);
}

