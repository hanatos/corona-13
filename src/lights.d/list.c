#include "lights.h"

// handles geometric light sources. nee to volumes and envmap is done in
// pathspace/nee instead.  this also has some (half-wired) code to start light
// paths on emitters at the bottom of this file. (TODO: wire volumes and sky)

typedef struct lights_t
{
  // light source importance sample list, adjustable
  float *prim_area;
  float *L;
  primid_t *primid;
  uint32_t num_alloced_prims;
  uint32_t num_prims;
  uint32_t inited;

  float p_geo, p_sky, p_vol;

  shader_so_t *vol;
}
lights_t;

lights_t *lights_init()
{
  lights_t *s = malloc(sizeof(*s));
  memset(s, 0, sizeof(*s));
  s->inited = 0;
  s->primid = 0;
  s->prim_area = 0;
  s->L = 0;
  s->num_alloced_prims = s->num_prims = 0;
  s->p_geo = s->p_sky = s->p_vol = 0.0f;
  return s;
}

void lights_cleanup(lights_t *s)
{
  free(s->L);
  free(s->prim_area);
  free(s->primid);
  free(s);
}

void lights_pdf_type(const path_t *p, const int v, mf_t p_egv[3])
{
  p_egv[0] = mf_set1(rt.lights->p_sky);
  p_egv[1] = mf_set1(rt.lights->p_geo);
  p_egv[2] = mf_set1(rt.lights->p_vol);
}

void lights_init_volume_light(shader_so_t *vol)
{
  rt.lights->vol = vol;
}

void lights_init_light(const uint32_t shapeid, const float L)
{
  uint32_t off = rt.lights->num_prims;
  rt.lights->num_prims += rt.prims->shape[shapeid].num_prims;
  if(rt.lights->num_alloced_prims < rt.lights->num_prims)
  {
    rt.lights->primid = realloc(rt.lights->primid, sizeof(primid_t)*rt.lights->num_prims);
    rt.lights->prim_area = realloc(rt.lights->prim_area, sizeof(float)*rt.lights->num_prims);
    rt.lights->L = realloc(rt.lights->L, sizeof(float)*rt.lights->num_prims);
    rt.lights->num_alloced_prims = rt.lights->num_prims;
  }
  for(int k=0;k<rt.prims->shape[shapeid].num_prims;k++)
  {
    rt.lights->primid[off+k] = rt.prims->shape[shapeid].primid[k];
    rt.lights->primid[off+k].shapeid = shapeid;
    rt.lights->prim_area[off+k] = prims_get_area(rt.prims, rt.lights->primid[off+k]) * L;
    rt.lights->L[off+k] = L;
  }
}

void lights_prepare_frame()
{
  if(rt.lights->inited) return;
  rt.lights->inited = 1;

  if(!rt.shader->skyshader.black) rt.lights->p_sky = 1.0f;
  if(rt.lights->num_prims) rt.lights->p_geo = 1.0f;
  if(rt.lights->vol) rt.lights->p_vol = 1.0f;
  float p_sum = rt.lights->p_sky + rt.lights->p_geo + rt.lights->p_vol;
  if(p_sum <= 0.0f) return;
  rt.lights->p_sky /= p_sum;
  rt.lights->p_geo /= p_sum;
  rt.lights->p_vol /= p_sum;
  if(rt.lights->num_prims == 0) return;

  float sum_area_L = 0.0f;
  for(int k=0;k<rt.lights->num_prims;k++)
    sum_area_L += rt.lights->prim_area[k];
  for(int k=0;k<rt.lights->num_prims;k++)
    rt.lights->L[k] /= sum_area_L;

  // calculate cdf
  for(int k=1;k<rt.lights->num_prims;k++)
    rt.lights->prim_area[k] += rt.lights->prim_area[k-1];
  //normalise
  for(int k=0;k<rt.lights->num_prims-1;k++)
    rt.lights->prim_area[k] /= rt.lights->prim_area[rt.lights->num_prims-1];
  rt.lights->prim_area[rt.lights->num_prims-1] = 1.0f;
}

mf_t lights_pdf_next_event(const path_t *p, int v)
{
  if(primid_invalid(p->v[v].hit.prim)) return mf_set1(0.0f);
  assert(v == 0 || v == p->length-1);
  uint32_t s = p->v[v].hit.prim.shapeid;
  // find shapeid in our list. this loop should be reasonably short (few of them in a scene)
  unsigned int min = 0, max = rt.lights->num_prims;
  unsigned int t = max/2;
  while (t != min)
  {
    if(rt.lights->primid[t-1].shapeid < s) min = t;
    else max = t;
    t = (min + max)/2;
  }
  // sanity check: have we been asked for the pdf of a non-emissive shape?
  if(rt.lights->primid[t].shapeid != s) return mf_set1(0.0f);
  // pdf = (L * prim_area)/sum_L_area * 1/prim_area
  //     = L/sum_L_area 
  // could also re-evaluate
  // shader->emission(&k, &em, 525.0f, self->data);
  // and only store normalisation weight
  return mf_set1(rt.lights->L[t]); // same for the same shape
}

mf_t _lights_sample_next_event(path_t *p, float r1, float r2, float r3)
{
  // sample light shader
  int v = p->length;

  // sample primitive
  unsigned int t = sample_cdf(rt.lights->prim_area, rt.lights->num_prims, r1);

  // sample point on triangle
  p->v[v].hit.prim = rt.lights->primid[t];
  prims_sample(rt.prims, rt.lights->primid[t], r2, r3, &p->v[v].hit, p->time);

  if(v)
  {
    // only init segment if we're not starting a light ray
    for(int k=0;k<3;k++)
      p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
    p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
    for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
  }
  // init emission and normals and such
  shader_prepare(p, v);
  
  p->v[v].pdf = mf_set1(rt.lights->L[t]);
  if(p->v[v].shading.roughness > 1.0f-1e-4f)
    p->v[v].material_modes = p->v[v].mode = s_emit | s_diffuse;
  else
    p->v[v].material_modes = p->v[v].mode = s_emit | s_glossy;
  return mf_div(p->v[v].shading.em, p->v[v].pdf);
}

// returns throughput into direction of segment, no visibility is checked here.
mf_t lights_sample_next_event(path_t *p)
{
  const int v = p->length;
  mf_t edf = _lights_sample_next_event(p, pointsampler(p, s_dim_nee_light2), pointsampler(p, s_dim_nee_x), pointsampler(p, s_dim_nee_y));
  if(p->v[v].shading.roughness > 1.0f-1e-4f)
    edf = mf_mul(edf, mf_set1(1.0f/M_PI));
  else
  {
    const float phongexp = 2.0f/(p->v[v].shading.roughness * p->v[v].shading.roughness) - 2.0f;
    edf = mf_mul(edf, mf_set1(powf(-dotproduct(p->v[v].hit.gn, p->e[v].omega), phongexp) * (phongexp + 2.0f)/(2.0f*M_PI)));
  }
  return edf;
}

// sample a light ray
mf_t lights_sample(path_t *p)
{
  assert(p->v[0].mode & s_emit);
  const mf_t throughput = _lights_sample_next_event(p, pointsampler(p, s_dim_lightsource), pointsampler(p, s_dim_light_x), pointsampler(p, s_dim_light_y));
  // now v[0] is a valid and prepare()d vertex

  // cos sample ray direction
  float x, y, z;
  const float phongexp = 2.0f/(p->v[0].shading.roughness * p->v[0].shading.roughness) - 2.0f;
  sample_cos_k(&x, &y, &z, phongexp+1, pointsampler(p, s_dim_edf_x), pointsampler(p, s_dim_edf_y));
  for(int k=0;k<3;k++)
    p->e[1].omega[k] = z*p->v[0].hit.n[k] + x*p->v[0].hit.a[k] + y*p->v[0].hit.b[k];
  normalise(p->e[1].omega);
  // p(y) = 1/A, p(omega) = cos^(k+1) (k+2)/2pi [dw], edf = cos^k(k+2)/2pi and one cos from lambert normalises it to 1 when integrating over it.
  // this is in projected solid angle measure [dwp], so a diffuse light source returns 1/M_PI
  mf_t pdf_omega = mf_set1(1.0f);
  if(p->v[0].shading.roughness > 1.0f-1e-4f)
  {
    pdf_omega = mf_set1(1.0f/M_PI);
    p->v[0].material_modes = p->v[0].mode = s_emit | s_diffuse;
  }
  else
  {
    pdf_omega = mf_set1(powf(z, phongexp) * (phongexp+2.0f)/(2.0f*M_PI));
    p->v[0].material_modes = p->v[0].mode = s_emit | s_glossy;
  }
  // v[0].pdf is set by sample_next_event
  p->v[1].pdf = mf_mul(pdf_omega, p->v[0].pdf);
  p->v[0].pdf = mf_set1(1.0f);
  p->v[0].rand_cnt = s_dim_num_lt_beg;
  p->v[1].rand_beg = p->v[0].rand_beg + p->v[0].rand_cnt;
  p->v[1].rand_cnt = 1; // only free path length

  return throughput;
}

mf_t lights_pdf(const path_t *p, int v)
{
  assert(v == 0 || v == p->length-1);
  // two cases: p->v[0] or p->v[p->length-1].
  const int e = (v==0) ? 1 : (p->length-1);
  const float *omega = p->e[e].omega;
  const hit_t *hit = &p->v[v].hit;
  const float z = dotproduct(hit->gn, omega);
  // no probability to shoot behind the ls
  if(p->v[0].mode & s_emit)
  {
    if(z <= 0.0f) return mf_set1(0.0f);
  }
  else if(p->v[0].mode & s_sensor)
  {
    if(z >= 0.0f) return mf_set1(0.0f);
  }

  const mf_t pdf_x = lights_pdf_next_event(p, v);

  const float phongexp = 2.0f/(p->v[v].shading.roughness * p->v[v].shading.roughness) - 2.0f;
  if(p->v[v].shading.roughness > 1.0f-1e-4f)
    return mf_mul(pdf_x, mf_set1(1.0f/M_PI));
  else
    return mf_mul(pdf_x, mf_set1(powf(fabsf(z), phongexp) * (phongexp+2.0f)/(2.0f*M_PI)));
}

// evaluate emitted radiance of given path vertex v,
// in direction towards the sensor.
mf_t lights_eval_vertex(path_t *path, int v)
{
  // early out for non-emissive vertices.
  // could have checked for mode & s_emit, but this is tighter:
  if(mf_all(mf_lte(path->v[v].shading.em, mf_set1(0.0f)))) return mf_set1(0.0f);
  mf_t edf = mf_set1(1.0f);
  // for path tracing, take edge incoming to light source
  const float *omega = path->e[v].omega;
  if(path->v[0].mode & s_emit) // for light tracing, use outgoing:
    omega = path->e[v+1].omega;
  // only geometric lights have phong edf with roughness, not envmaps nor volumetric emitters
  if(!primid_invalid(path->v[v].hit.prim))
  {
    // normal check for backfacing light sources:
    // we're only emitting in direction of geometric normals, but we need to distinguish between
    // path tracing and light tracing case.
    // this check is still valid when exiting a transmissive material that also
    // emits light to the outside (biolum on skin..):
    if((path->v[0].mode & s_emit) && dotproduct(path->v[v].hit.gn, omega) <= 0.0) return mf_set1(0.0f);
    else if((path->v[0].mode & s_sensor) && dotproduct(path->v[v].hit.gn, omega) >= 0.0) return mf_set1(0.0f);
    // eval phong edf
    if(path->v[v].shading.roughness > 1.0f-1e-4f)
    { // fast path for diffuse emission:
      edf = mf_set1(1.0f/M_PI);
    }
    else
    {
      const float phongexp = 2.0f/(path->v[v].shading.roughness * path->v[v].shading.roughness) - 2.0f;
      edf = mf_set1(powf(fabsf(dotproduct(path->v[v].hit.gn, omega)), phongexp) * (phongexp+2.0f)/(2.0f*M_PI));
    }
  }
  // envmap and volumetric lights only have isotropic emission
  return mf_mul(edf, path->v[v].shading.em);
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

