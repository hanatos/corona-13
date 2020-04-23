#pragma once

#include "shader.h"
#include "camera.h"
#include "pathspace/halfvec.h"

// this implements manifold next event estimation.
// it's a stubborn next event estimation technique, that doesn't give up easily.

// helper that returns 1 if a mnee connection is possible from vertex v
static inline int _mnee_connection_possible(const path_t *path, const int v)
{
  if(v < 1) return 0; // at least camera and first hit.
  if(v+1 >= PATHSPACE_MAX_VERTS) return 0; // can't add more vertices
  if(path->v[v].flags & s_environment) return 0; // already terminated (envmap vertex)
  // if there are no non-specular components, fail next event estimation.
  if(!(path->v[v].material_modes & (s_diffuse | s_glossy | s_volume)))
    return 0;
  // if(primid_invalid(path->v[v].hit.prim)) return 0; // don't connect from within media
  return 1;
}

static inline float _mnee_pdf_end(const path_t *path, int cullv, int v)
{
  if(path->length < 3) return 0.0f; // no next event for 2-vertex paths.
  assert(v == 0 || v == path->length-1);
  int e  = v ? v : v+1;    // edge leading up to the next event vertex, again v==0 is adjoint
  if(!_mnee_connection_possible(path, cullv)) return 0.0f;

  mf_t p_egv[3];
  lights_pdf_type(path, 0, p_egv);

  // light tracer connects to camera
  if(((path->v[0].mode & s_emit)>0) ^ (v==0))
    return view_cam_pdf_connect(path, v);
  else if((path->v[v].flags & s_environment) && (p_egv[0] > 0))
    return p_egv[0] * shader_sky_pdf_next_event(path, v) * path_G(path, e);
  else if(primid_invalid(path->v[v].hit.prim) && (path->v[v].mode & s_emit) && p_egv[2] > 0)
    return p_egv[2] * light_volume_pdf_nee(path, v);
  else if(!primid_invalid(path->v[v].hit.prim) && (path->v[v].mode & s_emit) && p_egv[1] > 0)
    return p_egv[1] * lights_pdf_next_event(path, v);
  return 0.0f;
}

static inline float _mnee_sample_end(path_t *p, int v)
{
  // constructing v[v] here:
  memset(p->v+v, 0, sizeof(vertex_t));
  memset(p->e+v, 0, sizeof(edge_t));
  p->v[v].pdf_mnee = -1.f;

  // instruct kelemen mlt to use new random numbers:
  p->v[v].rand_beg = p->v[v-1].rand_beg + p->v[v-1].rand_cnt;

  mf_t p_egv[3];
  lights_pdf_type(p, 0, p_egv);

  float edf = 0.0f;
  const float rand = pointsampler(p, s_dim_nee_light1);
  if(p->v[0].mode & s_emit)
  { // connect to the camera
    edf = view_cam_connect(p);
  }
  else if(rand < p_egv[0])
  { // connect to the envmap
    edf = shader_sky_sample_next_event(p);
    edf /= p_egv[0];
    // convert to vertex area measure
    const float G = path_G(p, v);
    p->v[v].pdf *= p_egv[0]*G;
  }
  else if(rand < p_egv[0] + p_egv[1])
  { // connect to geo lights
    edf = lights_sample_next_event(p);
    // compensate for envmap sampling probability
    p->v[v].pdf *= p_egv[1];
    edf /= p_egv[1];
  }
  else if(rand < p_egv[0] + p_egv[1] + p_egv[2])
  { // connect to volume lights
    edf = light_volume_sample_nee(p, v);
    p->v[v].pdf *= p_egv[2];
    edf /= p_egv[2];
  }
  else
  {
    fprintf(stderr, "[mnee] broken random number %g\n", rand);
  }
  return edf;
}

// create a nee postfix on the path.
static inline int _mnee_create_postfix(
    path_t *path,
    const int v,   // last vertex of fixed path, start next event estimation heer
    int vv)        // next event estimation vertex on the path. will be moved to the end.
{
  // store away the last vertex in case it's on the camera and propagation would destroy it:
  vertex_t nee_vertex = path->v[vv];
  // compute direction and maximum tracing distance:
  float dir[3] = {0.0};
  for(int k=0;k<3;k++) dir[k] = path->v[vv].hit.x[k] - path->v[v].hit.x[k];
  const float max_total_dist = sqrtf(dotproduct(dir, dir));
  for(int k=0;k<3;k++) dir[k] /= max_total_dist;
  float total_dist = 0.0f;

  // evaluate dummy bsdf at connecting vertex, to set v[v].mode needed by volume init for edge e[v+1]
  (void)shader_brdf(path, v);

  // keep direction and construct segments via path_propagate
  path->length = v+1;
  vv = v+1;
  for(;vv<PATHSPACE_MAX_VERTS;vv++)
  {
    // constructing v[vv] here:
    memset(path->v+vv, 0, sizeof(vertex_t));
    memset(path->e+vv, 0, sizeof(edge_t));
    // instruct kelemen mlt to use new random numbers:
    path->v[vv].rand_beg = path->v[vv-1].rand_beg + path->v[vv-1].rand_cnt;
    // init new outgoing direction
    for(int k=0;k<3;k++) path->e[vv].omega[k] = dir[k];

    // trick path_propagate() into not sampling a scattering point in volumes:
    // do that by pretending we're an mlt mutation and initing the next
    // vertex's primid to something not obviously invalid (will only be read once and then overwritten)
    path->v[vv].hit.prim = (primid_t){0, 0, 0, 0, 1};
    if(path_propagate(path, vv, s_propagate_mutate)) return -1;
    path->length++; // now v[vv] is a valid vertex
    assert(path->length == vv+1);
    total_dist += path->e[vv].dist;

    // TODO: optional optimisation: discard last vertex in case it doesn't have v[vv-1].material_modes & s_transmit
    // TODO: this should be a lot better in the presence of difficult occlusion, need to check occlusion only once at the end
    // TODO: need to remove material_modes check below!

    // passed the interesting point already:
    if(total_dist >= max_total_dist-rt.epsilon)
    { // found all interesting vertices, hooray :)
      if(path->v[vv].shading.em <= 0.0) return -1;

      // restore last vertex (in case of sensor we can't intersect, also want to restore pdf)
      path->v[vv] = nee_vertex;

      // implicitly: path->length = vv+1;
      // some bs ray tunneling fallback cases here
      if(!((path->v[vv].mode & s_sensor) || (path->v[vv].mode & s_emit)) ||
         !primid_eq(path->v[vv].hit.prim, nee_vertex.hit.prim))
          return -1;
      break;
    }
    // hit the environment map
    if(path->v[vv].flags & s_environment) return -1;

    // can only connect through transmissive interfaces
    if(!(path->v[vv].material_modes & s_transmit)) return -1;

    // set vertex interaction mode to something transmissive
    path->v[vv].mode = path->v[vv].material_modes & (s_diffuse|s_glossy|s_specular|s_transmit);
  }
  if(vv >= PATHSPACE_MAX_VERTS) return -1;
  assert(path->length == vv+1);
  return vv;
}

// helper function to fix a path according to fermat's law given a set of desired half vectors.
// returns dh_dx or 0.0 if failed
static inline int _mnee_fix_path(
    path_t *path,   // already filled path with additional T vertices and the next event vertex at the end (vv)
    const int v,    // last fixed vertex of prefix path      (start of chain)
    const int vv,   // next event vertex at very end of path (end of chain)
    const float *h, // newly sampled half vectors we're trying to converge to
    halfvec_stats_t *stats)
{
  if(vv-v <= 1) return 0;
  path->v[v].diffgeo.type = s_pinned_position;
  for(int k=v+1;k<vv;k++)
    path->v[k].diffgeo.type = s_free;

  // compute derivatives and half vectors:
  if(manifold_compute_tangents(path, v, vv) == 0.0) return 1;
  vertex_t nee_vertex = path->v[vv];

  const float ret = halfvec_to_worldspace(stats, path, h, v, vv);
  if(ret == 0.0f) return 1;

  // restore emission flags and pdf
  // but keep recomputed directionally dependent EDF
  path->v[vv].pdf = nee_vertex.pdf;
  path->v[vv].mode = nee_vertex.mode;
  path->v[vv].material_modes = nee_vertex.material_modes;
  return 0;
}

static inline float _mnee_stepsize(const float roughness)
{
  return roughness;
  // from looking at the beckmann formula, we would expect this:
  // return roughness * 1.0f/sqrtf(2.0f);
}

// TODO: replace by bsdf sampling and then reconstruct half vector
static inline void _mnee_sample_h(const float roughness, const float cos_in, const float r1, const float r2, float *h)
{
#if 1 // sample same lobe as phong rough dielectric would
  float fudge_factor = 1.2f - 0.2f * sqrtf(fabsf(cos_in));
  const float exponent_f = (2.0f/(fudge_factor * roughness * roughness) - 2.0f);
  float cosh;
  sample_cos_k(h, h+1, &cosh, exponent_f, r1, r2);
  // to plane/plane:
  h[0] /= cosh;
  h[1] /= cosh;
#else // gaussian in plane/plane
  const float stepsize = _mnee_stepsize(roughness);
  float g1, g2;
  sample_gaussian(r1, r2, &g1, &g2);
  h[0] = g1 * stepsize;
  h[1] = g2 * stepsize;
#endif
}

static inline float _mnee_pdf_h(const float roughness, const float cos_in, const float *h)
{
#if 1 // same lobe as phong rough dielectric would use
  float fudge_factor = 1.2f - 0.2f * sqrtf(fabsf(cos_in));
  const float exponent_f = (2.0f/(fudge_factor * roughness * roughness) - 2.0f);
  const float cosh = 1.0f/sqrtf(1.0f + h[0]*h[0] + h[1]*h[1]); // cosine of plane/plane input
  const float pdf_ph = powf(cosh, exponent_f) * (exponent_f + 1.0f)/(2.0f*M_PI);
  // convert to plane/plane
  return pdf_ph * cosh*cosh*cosh*cosh;
#else // gaussian in plane/plane
  const float stepsize = _mnee_stepsize(roughness);
  const float stdv0 = stepsize;
  const float stdv1 = stepsize;
  return 1.0f/(2.0f*M_PI * stdv0*stdv1) * exp(-0.5 * (h[0]*h[0]/(stdv0*stdv0) + h[1]*h[1]/(stdv1*stdv1)));
#endif
}

static inline double mnee_pdf(
    const path_t *path,
    const int mlength,
    halfvec_stats_t *stats);

// adds a number of vertices at the end, the last one being on the light source
// or camera pupil, depending on tracing direction of the path.
static inline int mnee_sample(path_t *path, halfvec_stats_t *stats)
{
  const int v = path->length-1; // last fixed vertex of old path

  if(!_mnee_connection_possible(path, v)) return 1;

  // sample point on other side, don't increment path length
  const float edf = _mnee_sample_end(path, v+1);
  if(edf <= 0.0f) goto fail;
  path->length = v+2; // include nee vertex in path length, or else we'll confuse path space functions.

  // we're sampling a direction on the environment light in dwp measure of the
  // current culling vertex (v).  we want to represent the pdf in dwp of the
  // last vertex in the chain later on, after we fixed the path.  so need to
  // track the ratio of cosines and correct it:
  double dwp_ratio = 1.0;
  if(path->v[v+1].flags & s_environment)
    dwp_ratio = 1.0/path_G(path, v+1);
  
  const int vv = _mnee_create_postfix(path, v, v+1);
  if(vv < 0) goto fail;
  
  // sample half vectors, put into this array (leave at 0 for specular constraints)
  // also need to divide out pdf for these further down in throughput computation! for non-specular materials
  // the respective jacobians are transparently switched on in halfvec_measurement().
  float h[2*PATHSPACE_MAX_VERTS] = {0.0f};
  float roughness[PATHSPACE_MAX_VERTS] = {0.0f/0.0f};
  float cos_in[PATHSPACE_MAX_VERTS];

  for(int k=v+1;k<vv;k++)
  {
    // pointsampler might decide random dimension based on path length, so
    // we need to re-adjust that:
    path->length = k;
    path->v[k].rand_beg = path->v[k-1].rand_beg + path->v[k-1].rand_cnt;
    // sample gaussians based on roughness
    if(path->v[k].mode & s_specular)
    {
      path->v[k].rand_cnt = 0;
    }
    else
    {
      roughness[k] = path->v[k].shading.roughness;
      cos_in[k] = dotproduct(path->e[k].omega, path->v[k].hit.n);
      // TODO: return pdf, too!
      _mnee_sample_h(roughness[k], cos_in[k], pointsampler(path, 0), pointsampler(path, 1), h + 2*(k-v));

      path->v[k].rand_cnt = 2;
    }
  }
  // add back last vertex, too
  path->length = vv+1;

  // will overwrite pdf of vertices due to projection and slightly offset half vectors
  if(_mnee_fix_path(path, v, vv, h, stats)) goto fail;

  float bsdf = shader_brdf(path, v);
  // compute the throughput of the glossy chain:
  // thr = measurement/pdf of the glossy chain (both in vertex area measure)
  // need to divide out half vector pdf (cancel with specular materials)
  // this is the half vector space measurement contribution (including determinant of transfer matrix).
  // the way we compute it it also contains volume transmittance.

  // need updated transfer matrix for measurement computation
  // if(raydifferentials_compute_rd_h_double(path, 0, v, vv)) goto fail;
  double dh_dx = raydifferentials_compute_rd_h(path, 0, v, vv);
  if(dh_dx == 0.0) goto fail;

  // bsdf called from within this sets reflect/transmit flags.
  double thr = halfvec_measurement(path, v, vv);

  double pdf_h = 1.0f;
  // now we have initialized half vectors, compute
  // the pdf for those (to precisely match the only approximate half vecs), also init
  // the jacobians for the conversion to vertex area pdf:
  for(int k=v+1;k<vv;k++)
  {
    // make sure we didn't accidentially flip sides:
    if(!(path->v[k].mode & s_transmit)) goto fail;
    if(!(path->v[k].mode & s_specular))
    {
      // TODO: either remove this and return pdf above, or do the same as in mnee_pdf!
      path->v[k].pdf = _mnee_pdf_h(roughness[k], cos_in[k], path->v[k].diffgeo.h);
      pdf_h *= path->v[k].pdf; // remember accumulated half vector space pdf
    }
    else path->v[k].pdf = 1.0f;
  }

  // convert inner vertices to vertex area measure:
  // multiply jacobian to random inner vertex (so it doesn't confuse the measurement below)
  // it is really not meaningful to assign individual vertex area pdf, these only go together.
  if(vv - v > 1) path->v[vv-1].pdf *= dh_dx;

  assert(vv==path->length-1);
  assert(path->v[vv].mode & s_emit);
  double throughput = path->v[v].throughput * bsdf * thr / (path->v[vv].pdf * pdf_h);
  // stupid clamping
#ifdef PTMNEE_BIAS
  throughput = MIN(throughput, 2.0);
#endif
  if(throughput <= 0.0) goto fail; // bsdf might be backfacing and 0 now.
  path->throughput = throughput;
  path->v[vv].total_throughput = path->throughput;
  if(path->v[vv].flags & s_environment) dwp_ratio *= path_G(path, vv);
  path->v[v].pdf_mnee = pdf_h * dh_dx * path->v[vv].pdf * dwp_ratio;
  assert(path->v[v].pdf_mnee > 0.0);
#if 0
  double check_pdf = mnee_pdf(path, v+1, stats);
  assert(fabs(path->v[v].pdf_mnee - check_pdf) < 2e-2*fmax(fmax(1.0, fabs(path->v[v].pdf_mnee)), fabs(check_pdf)));
#endif
  return 0;
fail:
  // roll back to old length
  path->length = v+1;
  return 1;
}

// returns the pdf of this path, assuming that the length of the path
// before doing mnee was the given length argument.
static inline double mnee_pdf(
    const path_t *path,
    const int mlength,
    halfvec_stats_t *stats)
{
  // pretty much the same code as the sample function above, only return pdf if path fixing worked.
  path_t tmp = *path;
  const int v = mlength-1;
  if(!_mnee_connection_possible(path, v)) return 0.0;

  // check pdf cache:
  if(path->v[v].pdf_mnee >= 0.0) return path->v[v].pdf_mnee;

  const float pdf_nee = _mnee_pdf_end(path, v, path->length-1);
  if(path->length - mlength == 1) return pdf_nee;

  // early out: check if it is a chain of transmissive events at all:
  for(int k=v+1;k<path->length-1;k++)
    if(!(path->v[k].material_modes & s_transmit)) return 0.0f;

  // could be optimized, only need to compute half vectors of the input path,
  // which might have been constructed by bsdf scattering:
  tmp.v[v].diffgeo.type = s_pinned_position;
  for(int i=v+1;i<tmp.length;i++)
    tmp.v[i].diffgeo.type = s_free;
  manifold_compute_tangents(&tmp, v, tmp.length-1);
  float h[2*PATHSPACE_MAX_VERTS] = {0.0f};
  for(int k=v+1;k<tmp.length-1;k++)
  {
    h[2*(k-v)+0] = tmp.v[k].diffgeo.h[0];
    h[2*(k-v)+1] = tmp.v[k].diffgeo.h[1];
  }

  // this will destroy the half vectors again:
  const int vv = _mnee_create_postfix(&tmp, v, path->length-1);
  if(vv < 0 || vv != path->length-1) return 0.0f;

  // store center path roughness (responsible for gaussian step size)
  float roughness[PATHSPACE_MAX_VERTS];
  float cos_in[PATHSPACE_MAX_VERTS];
  for(int k=v+1;k<vv;k++)
  {
    roughness[k-v] = tmp.v[k].shading.roughness;
    cos_in[k-v] = dotproduct(tmp.e[k].omega, tmp.v[k].hit.n);
  }

  // also mul half vector pdf, and convert to vertex area pdf.
  double pdf = 1.0;
  for(int k=v+1;k<vv;k++)
  {
    // get half vector pdf of all rough inner vertices
    if(!(path->v[k].mode & s_specular))
    {
      // TODO: compute pdf of sampling half vector of path->v[k] given
      // TODO: input direction of tmp->e[k]!
      pdf *= _mnee_pdf_h(roughness[k-v], cos_in[k-v], h + 2*(k-v));
    }
  }
  if(!(pdf > 0.0)) return 0.0f;
  if(_mnee_fix_path(&tmp, v, vv, h, stats)) return 0.0f;
#ifndef PTMNEE_BIAS
  // also need to check world space positions:
  for(int i=v+1;i<tmp.length-1;i++)
  {
    float dist[3];
    float quant = 1.0f;
    for(int k=0;k<3;k++) quant = fmaxf(fmaxf(fabsf(tmp.v[i].hit.x[k]), fabsf(path->v[i].hit.x[k])), quant);
    for(int k=0;k<3;k++) dist[k] = tmp.v[i].hit.x[k] - path->v[i].hit.x[k];
    const float len = sqrtf(dotproduct(dist, dist));
    if(len > quant * HALFVEC_REL_SPATIAL_EPS) return 0.0f;
  }
#endif

  if(vv - v > 1)
  { // need updated transfer matrix for measurement computation
    // if(raydifferentials_compute_rd_h_double(path, 0, v, vv)) goto fail;
    double dh_dx = raydifferentials_compute_rd_h(&tmp, 0, v, vv);
    if(dh_dx == 0.0) return 0.0;
    // convert inner vertices to vertex area measure:
    // this jacobian sometimes differs a little from the sampling one
    // unfortunately, due to slightly different tangent frames after
    // constructing the postfix.
    pdf *= dh_dx;
  }
  // TODO: in case of bdpt, we should store pdf_mnee here. need to invalidate all of them on every pop/extend though!
  return pdf_nee * pdf;
}

// pops the mnee created vertices from the end of the path, returning
// it to the state before calling mnee_next_event, restoring
// path->length to mlength.
static inline void mnee_pop(path_t *path, const int mlength)
{
  int v = mlength-1; // this is the new last valid vertex
  assert(mlength <= path->length);
  assert(mlength > 1);
  // TODO: this is probably not going to work with kelemen just like this (max rands/vert)
  // TODO: for all other vertices later on:
  for(int k=v+1;k<path->length;k++)
    path->v[v].rand_cnt += path->v[k].rand_cnt;
  // need to reset scatter mode, in case next vertex changes
  // from reflect to transmit. keep emission, though:
  path->v[v].mode &= s_emit;
  path->length = mlength;
  path->throughput = path->v[v].total_throughput; // reset throughput.
  path->v[v].pdf_mnee = -1.f;
}
