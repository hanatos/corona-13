#include "pathspace.h"
#include "points.h"
#include "spectrum.h"
#include "shader.h"
#include "prims.h"
#include "accel.h"
#include "view.h"
#include "pathspace/manifold.h"
#include "pathspace/tech.h"
#include "lights.h"
#include <string.h>

void path_init(path_t *path, uint64_t index, int camid)
{
  // init all required things to start a clean new path.
  path->lambda = mf_set1(0.0f);
  path->throughput = mf_set1(0.0f);
  path->length = 0;
  path->time = 0;
  path->temperature = 0;
  path->index = index;
  path->tangent_frame_scrambling = 0.0f;
  memset(&path->sensor, 0, sizeof(sensor_t));
  path->sensor.camid = camid;
  memset(path->v, 0, 2*sizeof(vertex_t));
  memset(path->e, 0, 2*sizeof(edge_t));
}

void path_copy(path_t *path, const path_t *other)
{
  path->length = other->length;
  path->lambda = other->lambda;
  path->time = other->time;
  path->temperature = other->temperature;
  path->throughput = other->throughput;
  path->index = other->index;
  path->tangent_frame_scrambling = other->tangent_frame_scrambling;
  path->sensor = other->sensor;
  for(int v=0;v<other->length;v++) path->v[v] = other->v[v];
  for(int e=0;e<other->length;e++) path->e[e] = other->e[e];
}

float path_lambert(const path_t *p, int v, const float *omega)
{
  if(!(p->v[v].mode & s_sensor) && primid_invalid(p->v[v].hit.prim)) return 1.0f; // no geo, no lambert (volumes + envmap)
  // fibers have lambert's law with sin of angle to tangent:
  if(p->v[v].material_modes & s_fiber)
  {
    const float cos_tw = dotproduct(p->v[v].diffgeo.dpdu, omega);
    return fmaxf(0.0f, 1.0f - cos_tw*cos_tw);
  }
  const float dot_sn = dotproduct(p->v[v].hit.n, omega);
  return fabsf(dot_sn);
}

// compute the regular geometric term for empty space transport
float path_G(const path_t *p, int e)
{
  if(p->v[e].flags & s_environment)
    return path_lambert(p, e-1, p->e[e].omega);
  if(p->v[e-1].flags & s_environment)
    return path_lambert(p, e, p->e[e].omega);
  return 
    path_lambert(p, e-1, p->e[e].omega)*
    path_lambert(p, e, p->e[e].omega)/
    (p->e[e].dist * p->e[e].dist);
}

// return the vertex number of the vertex with the interior that describes the
// volume of the requested edge e. if eta_ratio is set, assume transmit on v=e-1.
static inline int _path_edge_medium(const path_t *path, int e, int eta_ratio)
{
  int stack_v[PATHSPACE_MAX_VERTS] = {0};
  int sp = 1; // sp is number of stack entries, first is v[0] => global exterior medium
  for(int k=1;k<e;k++)
  {
    // assume transmit mode for last (k==v) vertex to find out ior of other side
    if(((k == e-1) && eta_ratio) || (path->v[k].mode & s_transmit))
    {
      if(!(path->v[k].flags & s_inside))
      { // entering a medium
        stack_v[sp++] = k; // push volume to list, order doesn't matter
      }
      else
      { // exiting a medium
        if(sp == 0) return -1.0f; // broken nesting, kill path
        for(int m=sp-1;m>=0;m--)
        {
          if(path->v[stack_v[m]].hit.prim.shapeid == path->v[k].hit.prim.shapeid)
          { // exited this medium, remove from list:
            stack_v[m] = stack_v[--sp];
            break;
          }
          if(m == 0) return -1.0f; //broken nesting.
        }
      }
    }
  }
  // find medium with highest priority (i.e. smallest shapeid, because invalid shapeid = (uint)-1)
  int result = 0;
  for(int i=1;i<sp;i++)
    if(((path->v[0].mode & s_emit) && stack_v[i]) ||  // for light tracing, always overwrite with anything not coming from v[0]. assumes light is in global exterior medium.
       (path->v[stack_v[i]].hit.prim.shapeid < path->v[result].hit.prim.shapeid))
      result = stack_v[i];
  return result;
}

mf_t path_eta_ratio(const path_t *path, int v)
{
  if(primid_invalid(path->v[v].hit.prim)) return mf_set1(1.0f); // ior doesn't change in medium
  const int mv = _path_edge_medium(path, v+1, 1);
  if(mv < 0) return mf_set1(-1.0f); // broken nesting
  // return eta ratio n1/n2 where n2 is on the other side of the medium.
  return mf_div(path->e[v].vol.ior, path->v[mv].interior.ior);
}

// init volume struct on edge path->e[v] by considering the transition at vertex path->v[v-1].
int path_edge_init_volume(path_t *path, int v)
{
  // volume boundary prepare() shaders always put the volume data into vertex.interior.
  // we now have to decide which volume to use for the next edge
  // if((path->length == 1) || (path->v[v-1].mode & (s_fiber|s_reflect|s_volume)))
  if(!(path->v[v-1].mode & s_transmit))
  {
    // stay in current volume if reflecting or volume scattering or fiber event
    assert(v > 0);
    path->e[v].vol = path->e[v-1].vol;
  } 
  else
  {
    const int mv = _path_edge_medium(path, v, 0);
    if(mv < 0) return 1;
    assert(mv >= 0 && mv < path->length);
    path->e[v].vol = path->v[mv].interior;
  }
  return 0;
}

void path_update_throughput(path_t* path, int v)
{
  // handle emissive edges: these need to be multiplied by the outgoing
  // throughput at v[v-1], but v.throughput stores the incoming throughput.
  // i.e. we want to do this before converting the outgoing throughput at v-1
  // to an incoming at v by multiplying transmittance / pdf below.
  path->throughput = mf_div(mf_mul(path->v[v].throughput, path->e[v].contribution), path->e[v].pdf);

  // get the rest of the measurement contribution which we did not sample (mostly 1/mu_t):
  path->v[v].throughput = mf_mul(path->v[v].throughput, mf_div(path->e[v].transmittance, path->e[v].pdf));

  // compute throughput of complete path. will evaluate to 0 if it's not complete.
  if((path->v[0].mode & s_sensor) && (path->v[v].mode & s_emit))
    path->throughput = mf_add(path->throughput, mf_mul(path->v[v].throughput, lights_eval_vertex(path, v)));
  else if((path->v[0].mode & s_emit) && (path->v[v].mode & s_sensor))
    path->throughput = mf_add(path->throughput, mf_mul(path->v[v].throughput, view_cam_eval(path))); // this technique is usually ignored though.
}

// sample new vertex. if length==0, sample new start point (from sensor or light)
int path_extend(path_t *path)
{
  // about to add vertex v == path->length should be < MAX_VERTS
  // no more memory, sorry guys.
  if(path->length >= PATHSPACE_MAX_VERTS) return 1;

  int v = path->length; // current vertex we're initing right now.

  // now initialising v[v] and e[v]. v==0 might have special stuff set from the outside, never
  // write over that. also, if called with len==0, we init two verts (done in path_init already)
  if(v)
  {
    memset(path->v+v, 0, sizeof(vertex_t));
    memset(path->e+v, 0, sizeof(edge_t));
    path->v[v-1].pdf_mnee = mf_set1(-1.f); // old vertex is now extended inner vertex
    path->v[v].pdf_mnee = mf_set1(-1.f);   // new vertex
  }

  // regular inner path vertex?
  if(path->length)
  {
    if(path->v[v-1].flags & s_environment) return 1;
    if(!mf_any(mf_gt(path->v[v-1].throughput, mf_set1(0.0f))))
    {
      path->v[v-1].throughput = mf_set1(0.0f);
      path->v[v-1].mode = s_absorb;
      return 1;
    }

    // 1) sample bsdf

    // remember our random number offset
    path->v[v].rand_beg = path->v[v-1].rand_beg + path->v[v-1].rand_cnt;
    path->v[v].pdf = mf_set1(1.0f); // everybody will just *= his sampling here.
    path->v[v].throughput = path->v[v-1].throughput;

    // sample the bsdf at vertex v-1. also inits e[v] with volume information.
    path->v[v].throughput = mf_mul(path->v[v].throughput, shader_sample(path));

    // this includes s_dim_russian_r, in case it is needed later
    // by means of calling path_russian_roulette(.).
    path->v[v].rand_cnt = s_dim_num_extend;
  }
  else
  {
    if(path->tangent_frame_scrambling == 0.0f)
      path->tangent_frame_scrambling = 0.1f + points_rand(rt.points, common_get_threadid())*(0.9f-0.1f);
    mf_t lambda_pdf;
    if(mf_all(mf_eq(path->lambda, mf_set1(0.0f))))
    { // if it's already set, don't sample a new one (for bdpt)
      // so do we stratify these now or not?
      float lf[MF_COUNT];
      for(int l=0;l<MF_COUNT;l++)
        lf[l] = fmodf(pointsampler(path, s_dim_lambda) + l/(float)MF_COUNT, 1.0f);
      path->lambda = spectrum_sample_lambda(mf_loadu(lf), &lambda_pdf);
      path->time = view_sample_time(path, pointsampler(path, s_dim_time)); // unit pdf, /= 1.0f to throughput.
    }

    // sample stereo camera, even for light paths
    path->sensor.camid = view_sample_camid(pointsampler(path, s_dim_camid));

    // we just started a new path, begin at random number 0.
    // path->v[0].rand_beg = 0; // this is done by path_init(). we don't overwrite it, since someone might have
                                // purposefully messed with it (kelemen mlt for bdpt)
    if(path->v[0].mode & s_emit)
      path->v[0].throughput = light_sampler_sample(path);
    else
      path->v[0].throughput = view_cam_sample(path);

    path->v[0].throughput = mf_div(path->v[0].throughput, mf_set1(view_pdf_camid(path->sensor.camid)));

    // sample start point always in global exterior medium:
    shader_exterior_medium(path);

    // we're actually creating two vertices at this point. the first is now done.
    path->length ++;
    // currently same for all the paths, so we ignore it
    // path->v[0].pdf *= lambda_pdf;
    // path->v[0].throughput /= lambda_pdf;
    path->v[1].throughput = path->v[0].throughput;
    path->v[0].tech = s_tech_extend;
    v++;
  }

  // propagate almost deterministically to next vertex
  // path->v[v]. takes care of volume distance sampling.
  if(mf_all(mf_lte(path->v[v].throughput, mf_set1(0.0f))) || path_propagate(path, v, s_propagate_sample))
  { // something went wrong (probably bsdf sampling).
    if(!(path->v[v-1].mode & s_emit)) // only set to absorption in case we didn't try to reflect off a light source.
      path->v[v-1].mode = s_absorb;   // the emit flag we'd like to keep.
    path->v[v].throughput = mf_set1(-0.0f);
    return 1; // return without incrementing path length.
  }

  // transform probability to on-surface probability at vertex v
  path->v[v].pdf = mf_mul(path->v[v].pdf, mf_set1(path_G(path, v)));
  path->length++; // now a valid vertex

  path_update_throughput(path, v);

  // store for potential path_pop()
  path->v[v].total_throughput = path->throughput;
  path->v[v].tech = s_tech_extend;
  return 0;
}

int path_russian_roulette(path_t *path, const float p_survival)
{
  assert(path->length); // need at least one vertex
  // perform russian roulette to terminate the path.
  // need to get random numbers for last vertex v:
  const int v = path->length - 1;
  path->length--;
  const float rr = pointsampler(path, s_dim_russian_r);
  path->length++;
  if(rr >= p_survival)
  { // die
    path->v[v].throughput = mf_mul(path->v[v].throughput, mf_set1(1.0f/(1.0f-p_survival)));
    path->v[v].pdf = mf_mul(path->v[v].pdf, mf_set1(1.0f-p_survival));
    return 1;
  }
  // survive
  path->v[v].throughput = mf_mul(path->v[v].throughput, mf_set1(1.0f/p_survival));
  path->v[v].pdf = mf_mul(path->v[v].pdf, mf_set1(p_survival));
  return 0;
}

void path_pop(path_t *path)
{
  assert(path->length > 2);
  int v = path->length-1;
  path->v[v-1].rand_cnt += path->v[v].rand_cnt;
  // need to reset scatter mode, in case next vertex changes
  // from reflect to transmit. keep emission, though:
  path->v[v-1].mode &= s_emit;
  path->length--;
  path->throughput = path->v[v-1].total_throughput; // reset throughput.

  // for mvnee: pop once more! vertex tech determines this:
  if(v && path->v[v-1].tech == s_tech_mvnee) path_pop(path);
}

// return 1 if p->v[v] is visible from p->v[v-1].
int path_visible(path_t *p, int v)
{
  ray_t ray;
  ray.time = p->time;
  float total_dist = prims_get_ray(&p->v[v-1].hit, &p->v[v].hit, &ray);
  // light tracer can't cull camera based on frame normal here.
  if(!(p->v[v].material_modes & s_volume) &&
     !(p->v[v].material_modes & s_fiber) &&
     !(p->v[0].mode & s_emit) && dotproduct(p->v[v].hit.gn, ray.dir) >= 0)
    return 0;
  hit_t hit = p->v[v-1].hit;
  // memset(&hit, 0, sizeof(hit));
  // hit.prim = p->v[v-1].hit.prim;
  vertex_t lightv = p->v[v]; // backup
  while(total_dist > 0.0f)
  {
    hit.dist = total_dist;
    hit.prim = INVALID_PRIMID;
    accel_intersect(rt.accel, &ray, &hit);
    if(hit.dist >= total_dist) break; // no intersection at all
    if(primid_invalid(hit.prim)) break; // gotta be the envmap here, volumes shouldn't happen
    if(primid_eq(hit.prim, lightv.hit.prim)) break; // reached the light
    total_dist -= hit.dist;
    // update hit and position for prepare()
    for(int k=0;k<3;k++)
      ray.pos[k] = hit.x[k] = ray.pos[k] + hit.dist * ray.dir[k];
    p->v[v].hit = hit;
    float prep = shader_prepare(p, v);
    if(prep >= 0.0) return 0;     // opaque object found
    prims_offset_ray(&hit, &ray); // only geo can be clip-mapped/non opaque
  }
  p->v[v] = lightv; // restore
  return 1;
}

md_t path_measurement_contribution_dwp(path_t *path, int s, int e)
{
#define MEASUREMENT_BEG s
#define MEASUREMENT_END e
#include "pathspace/measurement.h"
}

// get measurement contribution f in vertex area measure
md_t path_measurement_contribution_dx(path_t *path, int s, int e)
{
  // standard veach-style geometric terms:
#define MEASUREMENT_BEG s
#define MEASUREMENT_END e
#define PER_EDGE_JACOBIAN f = md_mul(f, md_set1(path_G(path, e)));
#include "pathspace/measurement.h"
}

// get pdf p as sampled, in product area measure. these numbers will be _very_ small.
md_t path_pdf(const path_t *path)
{
  // this works as long as we ignore the pdf of lambda sampling
  md_t pdf = md_set1(1.0);
  for(int k=0;k<path->length;k++)
    pdf = md_mul(pdf, mf_2d(path->v[k].pdf));
  return pdf;
}

// get throughput as sampled, X = f/p
mf_t path_throughput(const path_t *path)
{
  if(path->length < 2) return mf_set1(0.0f);
  if((path->v[0].mode & s_sensor) ||
     (path->v[path->length-1].mode & s_sensor))
    return path->throughput;
  return mf_set1(0.0f);
}

// return on-surface pdf of vertex v if it had been sampled via path_extend at path->length = v (from v[v-1])
mf_t path_pdf_extend(const path_t *path, int v)
{
  mf_t pdf = mf_set1(1.0f);
  if(v == 0)
    return mf_set1(1.0f); // constructed both at the same time.
  else if((v == 1) && (path->v[0].mode & s_emit))
  {
	  pdf = light_sampler_pdf_extend(path);
  }
  else if((v == 1) && (path->v[0].mode & s_sensor))
    pdf = view_cam_pdf(path, 0);
  else
    pdf = shader_pdf(path, v-1);
  // account for volume probabilities:
  mf_t vol_pdf = shader_vol_pdf(path, v);
  return mf_mul(mf_mul(vol_pdf, pdf), mf_set1(path_G(path, v)));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// extension of the reverse path from v+1
mf_t path_pdf_extend_adjoint(const path_t *path, int v)
{
  // last vertex is first adjoint one.
  // you can only create a single vertex the other way around
  // via next event estimation. extend will always produce at least two.
  if(v == path->length-1) return mf_set1(1.0f); // constructed both at the same time.
  
  mf_t pdf = mf_set1(1.0f);
  // special case for the first two vertices:
  if(v == path->length-2)
  {
    if(path->v[0].mode & s_emit)
    {
      // path is from light, so adjoint starts at camera
      // convert to area measure of v[l-2]
      pdf = view_cam_pdf(path, v+1);
    }
    else
    {
      pdf = light_sampler_pdf_extend_adjoint(path, v);
    }
  }
  else if ((v==0) && (path->v[0].mode & s_sensor))
  {
    // we don't intersect the camera aperture by chance with a light tracer
    return mf_set1(0.f);
  }
  else
  {
    // regular case
    pdf = shader_pdf_adj(path, v+1);
  }
  // account for volume probabilities:
  pdf = mf_mul(pdf, shader_vol_pdf_adjoint(path, v+1));
  // convert to vertex area of vertex after shooting the ray
  return mf_mul(pdf, mf_set1(path_G(path, v+1)));
}

void path_reverse(path_t *path, const path_t *input)
{
  if((input->length < 2) ||
    !((input->v[0].mode & s_sensor) || (input->v[input->length-1].mode & s_sensor)))
  {
    path_init(path, path->index, path->sensor.camid);
    return;
  }
  path->length = input->length;
  path->lambda = input->lambda;
  path->time   = input->time;
  path->temperature = input->temperature;
  path->throughput  = input->throughput;
  path->index  = input->index;
  path->sensor = input->sensor;

  const int n = input->length-1;
  for(int v=0;v<=n;v++)
    memcpy(path->v+v, input->v+n-v, sizeof(vertex_t));

  for(int e=1;e<=n;e++)
    memcpy(path->e+e, input->e+n+1-e, sizeof(edge_t));
  path->e[0] = path->e[1]; // copy over medium properties

  for(int k=0;k<=n;k++)
  { // flip direction:
    for(int i=0;i<3;i++)
      path->e[k+1].omega[i] = -path->e[k+1].omega[i];
    // take care of correct shading normal flipping and corresponding derivatives
    const int is_inside = path->v[k].flags & s_inside;
    const int should_be_inside = k && !primid_invalid(path->v[k].hit.prim) && !(path->v[k].mode & s_fiber) && (dotproduct(path->e[k].omega, path->v[k].hit.gn) > 0.0f);
    if(is_inside ^ should_be_inside)
    { // need to flip normal
      for(int i=0;i<3;i++) path->v[k].hit.n[i] *= -1.0f;
      for(int i=0;i<3;i++) path->v[k].diffgeo.dndu[i] *= -1.0f;
      for(int i=0;i<3;i++) path->v[k].diffgeo.dndv[i] *= -1.0f;
      path->v[k].flags ^= s_inside;
    }
    // keep track of technique:
    if     (path->v[k].tech == s_tech_extend)     path->v[k].tech = s_tech_extend_adj;
    else if(path->v[k].tech == s_tech_nee)        path->v[k].tech = s_tech_nee_adj;
    else if(path->v[k].tech == s_tech_extend_adj) path->v[k].tech = s_tech_extend;
    else if(path->v[k].tech == s_tech_nee_adj)    path->v[k].tech = s_tech_nee;
  }
  // seems to match, all good :)
  // double f1 = path_measurement_contribution_dwp(input);
  // double f2 = path_measurement_contribution_dwp(path);
  // fprintf(stderr, "measured %g %g\n", f1, f2);
  // assert(fabs(f1 - f2) < 1e-1f*fmax(f1, f2));
}

// connect two paths, extending path1 by a connection edge and the reverse of path2.
mf_t path_connect(path_t *path1, const path_t *path2)
{
  if(path1->length + path2->length > PATHSPACE_MAX_VERTS) return mf_set1(0.0f);
  // cut off paths ending on the envmap without connection:
  if(path1->v[path1->length-1].flags & s_environment) return mf_set1(0.0f);
  if(path2->v[path2->length-1].flags & s_environment) return mf_set1(0.0f);

  const int v = path1->length;

  // only connect same time and wavelength
  assert(path1->time == path2->time);
  assert(mf_all(mf_eq(path1->lambda, path2->lambda)));

  // init new segment before the loop below initing the inside flag:
  for(int i=0;i<3;i++)
    path1->e[v].omega[i] = path2->v[path2->length-1].hit.x[i] - path1->v[v-1].hit.x[i];
  path1->e[v].dist = sqrtf(dotproduct(path1->e[v].omega, path1->e[v].omega));
  for(int i=0;i<3;i++)
    path1->e[v].omega[i] /= path1->e[v].dist;
  path1->e[v].transmittance = mf_set1(0.0f);

  // append all vertices and segments of second path, in reverse order
  for(int k=0;k<path2->length;k++)
  {
    memcpy(path1->v+v+k, path2->v+path2->length-1-k, sizeof(vertex_t));
    memcpy(path1->e+v+k+1, path2->e+path2->length-1-k, sizeof(edge_t));
    // flip direction:
    for(int i=0;i<3;i++)
      path1->e[v+k+1].omega[i] = -path1->e[v+k+1].omega[i];
    // take care of correct shading normal flipping and corresponding derivatives
    const int is_inside = path1->v[v+k].flags & s_inside;
    const int should_be_inside = !primid_invalid(path1->v[v+k].hit.prim) && !(path1->v[v+k].mode & s_fiber) && (dotproduct(path1->e[v+k].omega, path1->v[v+k].hit.gn) > 0.0f);
    if(is_inside ^ should_be_inside)
    { // need to flip normal
      for(int i=0;i<3;i++) path1->v[v+k].hit.n[i] *= -1.0f;
      for(int i=0;i<3;i++) path1->v[v+k].diffgeo.dndu[i] *= -1.0f;
      for(int i=0;i<3;i++) path1->v[v+k].diffgeo.dndv[i] *= -1.0f;
      path1->v[v+k].flags ^= s_inside;
    }
    // keep track of technique:
    if     (path1->v[v+k].tech == s_tech_extend)     path1->v[v+k].tech = s_tech_extend_adj;
    else if(path1->v[v+k].tech == s_tech_nee)        path1->v[v+k].tech = s_tech_nee_adj;
    else if(path1->v[v+k].tech == s_tech_extend_adj) path1->v[v+k].tech = s_tech_extend;
    else if(path1->v[v+k].tech == s_tech_nee_adj)    path1->v[v+k].tech = s_tech_nee;
  }
  if(path1->v[0].mode & s_emit)
  {
    // take pixel position from eye path, if that's path2
    assert(!(path2->v[0].mode & s_emit));
    path1->sensor = path2->sensor;
  }

  // clear previous scatter flags:
  path1->v[v-1].mode &= s_emit;
  path1->v[v].mode &= s_emit;
  // eval bsdf to get modes for volume (T vs R)
  const mf_t bsdf1 = shader_brdf(path1, v-1);
  if(path_edge_init_volume(path1, v)) goto fail;
  const mf_t bsdf2 = shader_brdf(path1, v);
  if(mf_all(mf_lte(bsdf1, mf_set1(0.0f))) || mf_all(mf_lte(bsdf2, mf_set1(0.0f))) ||
    (path1->v[v-1].mode & s_specular) ||
    (path1->v[v].mode & s_specular))
    goto fail; // kills specular connections

  if(!path_visible(path1, v)) goto fail;

  mf_t vol = shader_vol_transmittance(path1, v);

  // fix cached eta ratios:
  for(int k=v;k<path1->length+path2->length;k++)
  {
    path1->v[k].diffgeo.eta = path_eta_ratio(path1, k);
    if(mf_all(mf_lt(path1->v[k].diffgeo.eta, mf_set1(0.0f)))) goto fail;
  }

  // evaluate throughput
  const mf_t connect = mf_mul(mf_mul(mf_mul(vol, mf_set1(path_G(path1, v))), bsdf1), bsdf2);
  const mf_t throughput = mf_mul(mf_mul(path1->v[v-1].throughput, connect), path2->v[path2->length-1].throughput);
  path1->length += path2->length; // mark edges and vertices as inited and in this path
  path1->throughput = throughput;
  return throughput;
fail:
  return mf_set1(0.0f);
}

// merge two paths in a biased way by joining the two end vertices into one.
// will update path1 and keep the joint vertex from path1.
// this is a deterministic, biased merge, if you're not doing this
// deterministically you should adjust the throughput by the pdf externally.
mf_t path_merge(path_t *path1, const path_t *path2)
{
  if(path1->length + path2->length > PATHSPACE_MAX_VERTS) return mf_set1(0.0f);
  // cut off paths ending on the envmap without connection:
  if(path1->v[path1->length-1].flags & s_environment) return mf_set1(0.0f);
  if(path2->v[path2->length-1].flags & s_environment) return mf_set1(0.0f);

  const int v = path1->length;

  for(int k=0;k<path2->length-1;k++)
  {
    memcpy(path1->v+v+k, path2->v+path2->length-2-k, sizeof(vertex_t));
    memcpy(path1->e+v+k+1, path2->e+path2->length-2-k, sizeof(edge_t));
    // flip direction:
    for(int i=0;i<3;i++)
      path1->e[v+k+1].omega[i] = -path1->e[v+k+1].omega[i];
    // take care of correct shading normal flipping and corresponding derivatives
    const int is_inside = path1->v[v+k].flags & s_inside;
    const int should_be_inside = !primid_invalid(path1->v[v+k].hit.prim) && !(path1->v[v+k].mode & s_fiber) && (dotproduct(path1->e[v+k].omega, path1->v[v+k].hit.gn) > 0.0f);
    if(is_inside ^ should_be_inside)
    { // need to flip normal
      for(int i=0;i<3;i++) path1->v[v+k].hit.n[i] *= -1.0f;
      for(int i=0;i<3;i++) path1->v[v+k].diffgeo.dndu[i] *= -1.0f;
      for(int i=0;i<3;i++) path1->v[v+k].diffgeo.dndv[i] *= -1.0f;
      path1->v[v+k].flags ^= s_inside;
    }
  }
  if(path1->v[0].mode & s_emit)
  {
    // take pixel position from eye path, if that's path2
    assert(!(path2->v[0].mode & s_emit));
    path1->sensor = path2->sensor;
  }

  // update segment (repoint to actual path)
  for(int i=0;i<3;i++)
    path1->e[v].omega[i] = path1->v[v].hit.x[i] - path1->v[v-1].hit.x[i];
  path1->e[v].dist = sqrtf(dotproduct(path1->e[v].omega, path1->e[v].omega));
  for(int i=0;i<3;i++)
    path1->e[v].omega[i] /= path1->e[v].dist;

  // clear previous scatter flags:
  path1->v[v-1].mode &= s_emit;
  // eval bsdf, get modes (T vs R)
  const mf_t bsdf = shader_brdf(path1, v-1);
  if(mf_all(mf_lte(bsdf, mf_set1(0.0f)))) goto fail; // kills specular connections

  // fix cached eta ratios:
  for(int k=v;k<path1->length+path2->length;k++)
    path1->v[v].diffgeo.eta = path_eta_ratio(path1, v);

  // evaluate throughput
  const mf_t throughput = mf_mul(mf_mul(path1->v[v-1].throughput, bsdf), path2->v[path2->length-1].throughput);
  path1->length += path2->length; // mark edges and vertices as inited and in this path
  path1->throughput = throughput;
  return throughput;
fail:
  return mf_set1(0.0f);
}

int path_project(path_t *p, int v, const path_propagation_mode_t mode)
{
  // if v[v-1] is the sensor, update pixel_{i,j}
  if(p->v[v-1].mode & s_sensor)
  {
    p->sensor.aperture_set = 1; // leave point on aperture as it was.
    p->v[v-1].throughput = view_cam_connect(p);
    if(mf_all(mf_lte(p->v[v-1].throughput, mf_set1(0.f)))) return 8;
  }
  else if((p->v[0].mode & s_emit) && (p->v[v].mode & s_sensor))
  {
    p->sensor.aperture_set = 1; // leave point on aperture as it was.
    p->v[v].throughput = view_cam_connect(p);
    if(mf_all(mf_lte(p->v[v].throughput, mf_set1(0.f)))) return 9;
  }

  // store projection position for reconstruction mode:
  const float oldx[3] = { p->v[v].hit.x[0], p->v[v].hit.x[1], p->v[v].hit.x[2] };
  const float eps = MAX(MAX(1.0f, fabsf(oldx[0])), MAX(fabsf(oldx[1]), fabsf(oldx[2])));
  const int old_envmap = p->v[v].flags & s_environment;

  // cast a ray from v[v-1] to v[v].hit.x and re-initialise the hitpoint at v
  // do that by just setting p->e[v].omega and calling path_propagate
  for(int k=0;k<3;k++) p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
  p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
  for(int k=0;k<3;k++) p->e[v].omega[k] /= p->e[v].dist;
  for(int k=0;k<3;k++) if(!(p->e[v].omega[k] == p->e[v].omega[k])) return 10; // dist near zero/volumes :(

  if((mode == s_propagate_reconstruct) && old_envmap)
    p->e[v].dist = FLT_MAX;
  if(mode == s_propagate_sample || (mode == s_propagate_mutate && !(p->v[v].mode & (s_sensor|s_volume))))
    p->e[v].dist = FLT_MAX;
  // if(v>1) // set mode correctly for updated direction.
    // shader_brdf(p, v-1); // we expect that to be done on the outside.

  const int err = path_propagate(p, v, mode);
  if(err) return err;

  if(mode == s_propagate_reconstruct || mode == s_propagate_mutate)
  { // test whether resulting vertex location is the one we asked for
    if(( (p->v[v].flags & s_environment) && !old_envmap) ||
       (!(p->v[v].flags & s_environment) &&  old_envmap))
      return 11; // pre and post envmap disagree
    if(mode == s_propagate_reconstruct)
    {
      if(p->v[v].flags & s_environment) return 0; // all good
      if(fabsf(p->v[v].hit.x[0] - oldx[0]) > 1e-3f*eps ||
         fabsf(p->v[v].hit.x[1] - oldx[1]) > 1e-3f*eps ||
         fabsf(p->v[v].hit.x[2] - oldx[2]) > 1e-3f*eps) return 12;
    }
  }
  return 0;
}

int path_propagate(path_t *path, int v, const path_propagation_mode_t mode)
{
  // take p->e[v].omega and shoot a ray from p->v[v-1] to init p->v[v].
  // take care of medium distance perturbation here if mutation_pdf != 0

  // 1) transition volume through interface using new direction.
  // if inconsistent fail the s_propagate_mutate call.
  if(path_edge_init_volume(path, v)) return 1;
  if((mode == s_propagate_mutate) && path->e[v].vol.shader == -1 && (path->v[v].mode & s_volume))
    return 13;

  // initialise to all bits set
  if(mode != s_propagate_mutate)
  { // mutations need to keep those as they were.
    // init to nothing, shader_sample() will set the scatter mode.
    path->v[v].mode = s_absorb;
    path->v[v].flags = s_none;
  }

  // distinguish between homo and hete media:
  // homo+scattering: first sample clip dist and then trace ray
  // hete or no scat: first trace ray and then sample dist up to there
  const int hete = shader_vol_hete(path, v);
  const int sample_vol_first = !hete && (mf(path->e[v].vol.mu_s, 0) > 0.0);

  // 2) cast ray/edge to get next vertex
  hit_t *hit = &path->v[v].hit;
  ray_t ray;

  // potentially clip distance at particle in volume:
  int vshader = -1;
  float clipdist = FLT_MAX;
  const float olddist = path->e[v].dist; // store for mutation
  if(mode == s_propagate_mutate)
  { // mutation and reconstruct cases keep exactly the distance to volume events as
    // sampled from the outside:
    if(primid_invalid(path->v[v].hit.prim))
      clipdist = path->e[v].dist;
  }
  else if(mode == s_propagate_reconstruct)
  { // always use the exact distance plus some terminator problem epsilon:
    const float eps = MAX(MAX(.5f, fabsf(path->v[v-1].hit.x[0])), MAX(fabsf(path->v[v-1].hit.x[1]), fabsf(path->v[v-1].hit.x[2])))*1e-3f;
    clipdist = path->e[v].dist + eps;
  }
  else if(sample_vol_first) // && s_propagate_sample
  { // for heterogeneous it is more efficient to trace the ray first.
    path->e[v].dist = FLT_MAX;
    clipdist = shader_vol_sample(path, v);
  }
  // assume we will scatter in the medium:
  if(clipdist < FLT_MAX) vshader = path->e[v].vol.shader;
  hit->prim = INVALID_PRIMID;
  hit->dist = clipdist;
  hit->shader = vshader;

  // cast ray:

  for(int k=0;k<3;k++) ray.pos[k] = path->v[v-1].hit.x[k];
  for(int k=0;k<3;k++) ray.dir[k] = path->e[v].omega[k];
  ray.time = path->time;
  ray.ignore = INVALID_PRIMID;
  ray.min_dist = 0.0;
  // offset only if on geometry:
  if(!primid_invalid(path->v[v-1].hit.prim))
    prims_offset_ray(&path->v[v-1].hit, &ray);

  accel_intersect(rt.accel, &ray, hit);
  float total_dist = hit->dist;
  path->e[v].dist = total_dist;

  // light tracer exited to envmap:
  if(hit->dist == FLT_MAX && (path->v[0].mode & s_emit) && !(path->v[v].mode & s_sensor))
    return 3;

  // reconstruction failed to find a vertex up until clip distance, but geometry was requested
  if(mode == s_propagate_reconstruct && hit->dist == clipdist && clipdist < FLT_MAX && path->e[v].vol.shader == -1)
    return 4;

  // 3) call prepare shader on new vertex, maybe have to trace more rays
  //    for shell primitives which are partially see-through
  if(hit->dist < FLT_MAX)
  {
    // update hit position:
    for(int k=0;k<3;k++)
      path->v[v].hit.x[k] = ray.pos[k] + hit->dist * ray.dir[k];
    // every new sample point not on the env map needs to call prepare() shaders:

    // TODO: this is only needed if clip maps are active (else we won't collect any vertices in the clipped distance)
    // TODO: maybe need to reset shading info here first (overwriting homo medium?)
    // TODO: if sample_vol_first store vertex (esp. shading + interior)
    float prep = shader_prepare(path, v);
    // but only geometric intersections can make ray tracing restart:
    while(hit->dist < FLT_MAX && prep < 0.0f && !primid_invalid(hit->prim))
    {
      for(int k=0;k<3;k++) ray.pos[k] += ray.dir[k] * hit->dist;
      clipdist -= hit->dist;
      hit->dist = clipdist;
      hit->shader = vshader;
      ray.ignore = hit->prim;
      if(!primid_invalid(hit->prim)) prims_offset_ray(hit, &ray);
      hit->prim = INVALID_PRIMID;
      accel_intersect(rt.accel, &ray, hit);
      total_dist += hit->dist;
      path->e[v].dist = total_dist;
      for(int k=0;k<3;k++)
        path->v[v].hit.x[k] = ray.pos[k] + hit->dist * ray.dir[k];
      if(hit->dist == FLT_MAX) break; // envmap doesn't prepare shaders
      prep = shader_prepare(path, v);
    }
    // TODO: maybe need to restore shading info from homo medium here, if clip maps
    //       didn't result in shorter distance!
    // if sample_vol_first and dist
  }

  // detect self-intersection, try not to detect sphere/cylinder TT paths.
  // should probably have a callback for convex objects instead of the magic > 2 test.
  if((path->v[v].hit.prim.vcnt > 2 || path->e[v].dist < 1e-4f) &&
      !primid_invalid(path->v[v].hit.prim) && primid_eq(path->v[v].hit.prim, path->v[v-1].hit.prim))
  {
#if 0
    fprintf(stderr, "detected self intersection at verts %d/%d\n", v-1, v);
    path_print(path, stderr);
    fprintf(stderr, "dot %f\n", dotproduct(path->v[v-1].hit.gn, path->e[v].omega));
    fprintf(stderr, "dist %f\n", path->e[v].dist);
    fprintf(stderr, "shape %d vcnt %d\n", path->v[v-1].hit.prim.shapeid, path->v[v-1].hit.prim.vcnt);
#endif
    return 5;
  }

  if(mode == s_propagate_sample)
  { // only evaluate transmittance/sample distance if we're not in mutation mode:
    if(!sample_vol_first)
    { // shader_vol_sample clipped to e[v].dist and return pdf, transmittance, emission:
      clipdist = shader_vol_sample(path, v);
      if(clipdist < path->e[v].dist)
      {
        hit->prim = INVALID_PRIMID;
        hit->shader = vshader;
        hit->dist = clipdist;
        path->e[v].dist = clipdist;
        for(int k=0;k<3;k++)
          path->v[v].hit.x[k] = path->v[v-1].hit.x[k] + hit->dist * path->e[v].omega[k];
        shader_prepare(path, v);
      }
    }
    else
    { // homo needs to update values based on sampled distance clipped to geometry
      if(path->e[v].vol.shader >= 0 && !(path->v[v].material_modes & s_volume))
        path->e[v].pdf = path->e[v].transmittance = shader_vol_transmittance(path, v);
    }
  }
  else if(mode == s_propagate_reconstruct)
  {
    // initialise volume attenuation:
    if(path->e[v].vol.shader >= 0)
    {
      for(int k=0;k<3;k++) // reconstruct precise location since we added an epsilon above
        path->v[v].hit.x[k] = path->v[v-1].hit.x[k] + olddist * path->e[v].omega[k];
      path->e[v].pdf = mf_set1(1.0f);
      path->e[v].transmittance = shader_vol_transmittance(path, v);
    }
  }

  if(path->e[v].dist >= FLT_MAX)
  { // envmap hit. this doesn't get shader_prepare() calls, in particular it doesn't have an inited manifold system
    path->v[v].flags |= s_environment;
    // init envmap emission
    shader_prepare(path, v);

    // just move out of the box a little, so we can compute tangent space half
    // vector derivatives with it:
    const float *aabb = accel_aabb(rt.accel);
    const float far = 2.0f*MAX(aabb[5] - aabb[2],
                           MAX(aabb[4] - aabb[1],
                               aabb[3] - aabb[0]));
    for(int k=0;k<3;k++)
      path->v[v].hit.x[k] = path->v[v-1].hit.x[k] + far * path->e[v].omega[k];
  }

  // set emission flags on all sorts of light sources (geo, vol, segment, and envmap)
  if(mf_any(mf_gt(path->e[v].contribution, mf_set1(0.0f))) ||
      (mf_any(mf_gt(path->v[v].shading.em, mf_set1(0.0f))) && !(path->v[v].flags & s_inside)))
    path->v[v].material_modes = path->v[v].mode = s_emit;

  if(mode == s_propagate_mutate)
  { // mutation mode
    if((path->v[0].mode & s_emit) && (path->v[v].mode & s_sensor))
    { // if connecting to the camera and intersecting something before we're quite there:
      if(hit->dist < olddist) return 6;
    }
    path->e[v].transmittance = shader_vol_transmittance(path, v);
  }

  // include distance sampling in projected solid angle pdf on vertex (will be
  // converted to on-surface on the outside, for instance in path_extend)
  if(mode == s_propagate_sample)
    path->v[v].pdf = mf_mul(path->v[v].pdf, path->e[v].pdf);

  return 0;
}

// printing related stuff
static inline char path_vertex_scatter_mode(const path_t *path, int v)
{
  if(path->v[v].mode & s_specular) return 'S';
  if(path->v[v].mode & s_diffuse)  return 'D';
  if(path->v[v].mode & s_glossy)   return 'G';
  if(path->v[v].mode & s_volume)   return 'V';
  return ' ';
}

static inline char path_vertex_kind(const path_t *path, int v)
{
  if(path->v[v].mode & s_reflect)  return 'R';
  if(path->v[v].mode & s_transmit) return 'T';
  if(path->v[v].mode & s_volume)   return 'o';
  if(path->v[v].mode & s_fiber)    return 'F';
  if(path->v[v].flags & s_environment) return 'I';
  if(path->v[v].mode & s_emit)     return 'L';
  if(path->v[v].mode & s_sensor)   return 'E';
  if(path->v[v].mode == s_absorb)  return 'X';
  return '?';
}

void path_print(const path_t *path, FILE *f)
{
  fprintf(f, "path at time %f wavelength %f pixel (%f, %f):\n", path->time, mf(path->lambda, 0), path->sensor.pixel_i, path->sensor.pixel_j);
  for(int k=0;k<path->length;k++)
  {
    if(path->e[k].vol.shader >= 0)
      fprintf(f, " |  tau %f em %f pdf %f dist %f ior %f\n", mf(path->e[k].transmittance, 0), mf(path->e[k].contribution, 0), mf(path->e[k].pdf, 0), path->e[k].dist, mf(path->e[k].vol.ior, 0));
    else
      fprintf(f, "    dist %f ior %f\n", path->e[k].dist, mf(path->e[k].vol.ior, 0));
    fprintf(f, " %c%c [%d] throughput %f eta %f\n",
        //(path->e[k].shader >= 0 || path->e[k+1].shader >= 0) ?
        path_vertex_kind(path, k),
        path_vertex_scatter_mode(path, k),
        k, mf(path->v[k].throughput, 0),
        mf(path_eta_ratio(path, k), 0));
  }
}

void path_print_manifold_matrix(const path_t *path, FILE *f)
{
  for(int k=0;k<path->length;k++)
  { // draw box outlines
    const int s = (path->v[k].mode & s_volume ? 3 : 2);
    for(int kk=0;kk<s;kk++)
    {
      if(k==0 && kk==0) fprintf(f, "\u250c");
      else              fprintf(f, "\u2500");
      for(int kkk=0;kkk<5;kkk++) fprintf(f, "\u2500");
      if(k==path->length-1 && kk==s-1) fprintf(f, "\u2510");
      else if(kk==s-1)                 fprintf(f, "\u252c");
      else                             fprintf(f, "\u2500");
    }
  }
  fprintf(f, "\n");
  for(int j=1;j<path->length-1;j++)
  {
    int rsc = path->v[j].mode & s_volume ? 3 : 2;
    for(int row=0;row<rsc;row++)
    {
      for(int k=1;k<j;k++)
        for(int kk=0;kk<(path->v[k-1].mode & s_volume ? 3 : 2);kk++)
          fprintf(f, "       ");
      int sp = path->v[j-1].mode & s_volume ? 3 : 2;
      int sc = path->v[j  ].mode & s_volume ? 3 : 2;
      int sn = path->v[j+1].mode & s_volume ? 3 : 2;
      for(int i=0;i<sp;i++)
        fprintf(f, "% 4.2f  ", path->v[j].diffgeo.a[sp*row+i]);
      for(int i=0;i<sc;i++)
        fprintf(f, "% 4.2f  ", path->v[j].diffgeo.b[sc*row+i]);
      for(int i=0;i<sn;i++)
        fprintf(f, "% 4.2f  ", path->v[j].diffgeo.c[sn*row+i]);
      fprintf(f, "\n");
    }
    if(j<path->length-2)
    {
      for(int k=0;k<path->length;k++)
      { // draw box outlines in between blocks
        const int s = (path->v[k].mode & s_volume ? 3 : 2);
        for(int kk=0;kk<s;kk++)
        {
          if(k==0 && kk==0) fprintf(f, "\u251c");
          else              fprintf(f, "\u2500");
          for(int kkk=0;kkk<5;kkk++) fprintf(f, "\u2500");
          if(k==path->length-1 && kk==s-1) fprintf(f, "\u2524");
          else if(kk==s-1)                 fprintf(f, "\u253c");
          else                             fprintf(f, "\u2500");
        }
      }
      fprintf(f, "\n");
    }
  }
  for(int k=0;k<path->length;k++)
  { // draw final box outline
    const int s = (path->v[k].mode & s_volume ? 3 : 2);
    for(int kk=0;kk<s;kk++)
    {
      if(k==0 && kk==0) fprintf(f, "\u2514");
      else              fprintf(f, "\u2500");
      for(int kkk=0;kkk<5;kkk++) fprintf(f, "\u2500");
      if(k==path->length-1 && kk==s-1) fprintf(f, "\u2518");
      else if(kk==s-1)                 fprintf(f, "\u2534");
      else                             fprintf(f, "\u2500");
    }
  }
  fprintf(f, "\n");
}

