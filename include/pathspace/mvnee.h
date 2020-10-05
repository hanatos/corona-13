#pragma once
#include "pathspace.h"
#include "pathspace/tech.h"
#include "shader.h"
#include "lights.h"
#include "view.h"

#include <math.h>

static inline double t_pdf(double t, double cosh)
{
  double theta = acosf(CLAMP(cosh, -1.0, 1.0));
  double sth = sin(theta);
  return sth/(theta * sqrt(1.0 - (1.0-2.0*t)*(1.0-2.0*t) * sth*sth));

#if 0
// FIXME: this is -nan at 0 and 1, though it seems it goes to a real value relatively flat!
  const double sinh = sqrt(1.0 - cosh*cosh);
  const double sinh2 = sinh*sinh;
  const double sinh3 = sinh2*sinh;
  const double sinh4 = sinh2*sinh2;
  const double t2 = t*t;
  return
  (cosh*sinh*
   cosh // sqrt(1.0-sinh2)
             * sqrt((-4.0*sinh2*t2)+4.0*sinh2*t-sinh2+1.0)
 +2.0*cosh*sinh3*t2-2*cosh*sinh3*t+cosh*sinh3-cosh*sinh)
 /(sqrt(1.0-sinh2)*sqrt((-4.0*sinh2*t2)+4.0*sinh2*t-sinh2+1.0)
                   *(2.0*sinh2*t2-2.0*sinh2*t+sinh2-1.0)
  +(4.0*sinh4-4*sinh2)*t2+(4*sinh2-4*sinh4)*t+sinh4-2*sinh2
  +1.0);
#endif
}

static inline float t_rand(uint64_t *state)
{
  uint64_t s1 = state[0];
  uint64_t s0 = state[1];
  state[0] = s0;
  s1 ^= s1 << 23;
  s1 ^= s1 >> 17;
  s1 ^= s0;
  s1 ^= s0 >> 26;
  state[1] = s1;
  // return (state0 + state1) / ((double)((uint64_t)-1) + 1.0);
  uint32_t v = 0x3f800000 | ((state[0]+state[1])>>41); // faster than double version.
  return (*(float*)&v) - 1.0f;
}

static inline double t_sample(double cosh, float r0)
{
  if(cosh > 0.9999) return 1.0;
  double theta = acosf(CLAMP(cosh, -1.0, 1.0));
  return 0.5*(1.0 + sin(theta *(2.0*r0-1.0))/sin(theta));
#if 0
#if 0
  for(int k=0;k<100;k++)
  {
    fprintf(stderr, "%g %g\n", k/990.0, t_pdf(k/990.0, cosh));
  }
  exit(0);
#endif
  // most stupid rejection sampling i could think of.
  // the function is max at t=0 and t=1, so we use this as bound:
  double max = t_pdf(1e-6, cosh);
  uint64_t state[2] = {(666ul<<32)*r0, 1337ul<<32};
  for(int i=0;i<100;i++) // limit max tries
  {
    float x = t_rand(state);
    float y = t_rand(state);
    float pdf = t_pdf(x, cosh);
    if(y*max < pdf) return x;
  }
  return -1;
#endif
}

static inline int
mvnee_possible(const path_t *p, const int v)
{
#if 0 // TODO: implement this. i think this is broken
  // if there are no non-specular components, fail next event estimation.
  if(!(p->v[v].material_modes & (s_diffuse | s_glossy)))
    return 0;
  // need to be inside some volume
  if(!primid_invalid(p->v[v].hit.prim)) return 0; // XXX but v is on the whale!
#endif
  return 1;
}

// this returns the product vertex area measure pdf of sampling
// both end vertices with mvnee.
static inline mf_t
mvnee_pdf(const path_t* p, int v)
{
  if(p->length < 3) return mf_set1(0.0f); // no next event for 2-vertex paths.
  assert(v == 0 || v == p->length-1);
  int v0 = v ? v-2 : v+2;  // vertex where next event was sampled from (v==0 means adjoint pdf)
  int v1 = v ? v-1 : v+1;  // middle vertex created by mvnee
  int e0 = v ? v-1 : 2;
  int e1 = v ? v   : 1;
  if(!mvnee_possible(p, v0)) return mf_set1(0.0f);

  mf_t p_egv[3];
  lights_pdf_type(p, v0, p_egv);

  mf_t res = mf_set1(0.0f);
  float cos_theta = dotproduct(p->e[e0].omega, p->e[e1].omega);
  if(cos_theta <= 0.0f) return mf_set1(0.0f);

  // light tracer connects to camera
  if(((p->v[0].mode & s_emit)>0) ^ (v==0))
    res = view_cam_pdf_connect(p, v);
  else if(p->v[v].flags & s_environment)
    res = mf_mul(p_egv[0], shader_sky_pdf_next_event(p, v));
  else if(primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf_any(mf_gt(p->v[v].shading.em, mf_set1(0.0f))) && mf(p_egv[2], 0) > 0)
    res = mf_mul(p_egv[2], light_volume_pdf_nee(p, v));
  else if(!primid_invalid(p->v[v].hit.prim) && (p->v[v].mode & s_emit) && mf(p_egv[1], 0) > 0)
    res = mf_mul(p_egv[1], lights_pdf_next_event(p, v));
  else return mf_set1(0.0f);

  // eval pdf as product of v and v1:
  float d[3] = {
    p->v[v].hit.x[0] - p->v[v0].hit.x[0],
    p->v[v].hit.x[1] - p->v[v0].hit.x[1],
    p->v[v].hit.x[2] - p->v[v0].hit.x[2]};
  float d1[3] = {
    p->v[v1].hit.x[0] - p->v[v0].hit.x[0],
    p->v[v1].hit.x[1] - p->v[v0].hit.x[1],
    p->v[v1].hit.x[2] - p->v[v0].hit.x[2]};
  float s = sqrtf(dotproduct(d, d));
  float t = dotproduct(d1, d) / (s * s);
  float hg_pdf = sample_eval_hg_fwd(p->v[v1].interior.mean_cos, p->e[e0].omega, p->e[e1].omega);   // should probably eval the phase function via shader_brdf()
  return mf_mul(res, mf_set1(2*M_PI/s * t_pdf(t, cos_theta) * hg_pdf));
}

// return on-surface pdf of vertex v if it had been sampled the other way around via
// next event estimation of the reverse path from v+1
static inline mf_t
mvnee_pdf_adjoint(const path_t *path, int v)
{
  // luckily these guys know already about v==0 or v==length-1
  return mvnee_pdf(path, v);
}

// sample next event. returns != 0 on failure
static inline int
mvnee_sample(path_t *p)
{
  assert(p->length >= 2); // at least camera and first hit.
  if(p->v[p->length-1].flags & s_environment) return 1;
  if(p->length >= PATHSPACE_MAX_VERTS-1) return 1; // need to append two new vertices
  const int v = p->length; // constructing new vertex here

  if(!mvnee_possible(p, v-1)) goto fail;

  // sample vertex v as endpoint (lights or camera)
  mf_t edf = mf_set1(0.0f);
  mf_t p_egv[3];
  lights_pdf_type(p, v-1, p_egv);

  // constructing v[v] here:
  memset(p->v+v, 0, sizeof(vertex_t));
  memset(p->e+v, 0, sizeof(edge_t));
  p->v[v].pdf_mnee = mf_set1(-1.f);

  // instruct kelemen mlt to use new random numbers:
  p->v[v].rand_beg = p->v[v-1].rand_beg + p->v[v-1].rand_cnt;

  // set technique to next event estimation
  p->v[v].tech = s_tech_mvnee;

  const float rand = pointsampler(p, s_dim_nee_light1);
  if(p->v[0].mode & s_emit)
  { // connect to the camera
    edf = view_cam_connect(p);
  }
  else if(rand < mf(p_egv[0],0))
  { // connect to the envmap
    edf = shader_sky_sample_next_event(p);
    edf = mf_div(edf, mf_set1(mf(p_egv[0], 0)));
    // pdf is in solid angle
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[0]);
  }
  else if(rand < mf(p_egv[0],0)+mf(p_egv[1],0))
  { // connect to geo lights
    edf = lights_sample_next_event(p);
    // compensate for envmap sampling probability
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[1]);
    edf = mf_div(edf, mf_set1(mf(p_egv[1], 0)));
  }
  else if(mf(p_egv[2],0) > 0)
  { // connect to volume lights
    edf = light_volume_sample_nee(p, v);
    if(!mf_any(mf_gt(edf, mf_set1(.0f)))) goto fail;
    p->v[v].pdf = mf_mul(p->v[v].pdf, p_egv[2]);
    edf = mf_div(edf, mf_set1(mf(p_egv[2],0)));
    // init segment
    for(int k=0;k<3;k++)
      p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
    p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
    for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
  }

  // ask edf and bsdf for their consent:
  if(!mf_any(mf_gt(edf, mf_set1(0.0f)))) goto fail;
  // compute brdf and throughput:
  mf_t bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1
  if(!mf_any(mf_gt(bsdf, mf_set1(0.0f)))) goto fail; // check for specular materials.

  // determine side of surface and volume from that (brdf sets mode)
  if(path_edge_init_volume(p, v)) goto fail;

  // move vertex v to v+1, free up space for v[v], our in-between vertex.
  p->e[v+1] = p->e[v];
  p->v[v+1] = p->v[v];
  // get volume properties:
  if(p->e[v].vol.shader == -1) goto fail; // no volume no mvnee.
  const float g = p->e[v].vol.mean_cos;
  // fprintf(stderr, "mean cos %g\n", g);
  // sample new vertex v on edge between the two:
  p->length = v+1; // instruct pointsampler to get new dimensions
  float hg_r1  = pointsampler(p, s_dim_nee_x); // used for cosh
  float hg_r2  = pointsampler(p, s_dim_nee_y); // usef for isotropic phi
  // sample hg with g
  float hg_pdf = 0.0f, omega[3];
  sample_hg_fwd(g, hg_r1, hg_r2, omega, &hg_pdf);
#if 0
  // can't do backscattering
  if(omega[0] < 0.0)
  {
    // rejection sample until we got it:
    uint64_t state[2] = {p->index, 1337*hg_r1};
    for(int k=0;k<100;k++)
    {
      hg_r1 = t_rand(state);
      sample_hg(g, hg_r1, hg_r2, omega, &hg_pdf);
      if(omega[0] >= 0.0f) break;
    }
  }
#endif

  // the cosine from phase function sampling defines the angle between e[v].w and e[v+1].w.
  const float dist_r = pointsampler(p, s_dim_nee_light1);
  const float t = dist_r; // XXX t_sample(omega[0], dist_r); // sampling is normalised to [0,1]

  // angle around axis of symmetry
  const float phi = 2.0f*M_PI*hg_r2;
  const float sin_phi= sinf(phi);
  const float cos_phi= cosf(phi);
  const float sinh2 = 1.0-omega[0]*omega[0];
  // shit, i think this may be subject to catastrophic cancellation
  // const float r = p->e[v].dist * (sqrt(1.0/(sinh2*4) - (0.5-t)*(0.5-t)) - sqrt(1.0/(sinh2*4)-1./4.));
  const double r = sqrt(1.0/(4.0*sinh2) - (0.5-t)*(0.5-t)) - sqrt(1.0/(4.0*sinh2) - 0.25);
  const double s = p->e[v].dist;

#if 0
  fprintf(stdout, "%g %g %g\n", t, r/p->e[v].dist, omega[0]);
  // gnuplot>
  // set size ratio -1
  // set yrange [0:0.4]
  // plot 'dat' u 1:2:3 w p ps 3 pt 7 lc palette
#endif

  // initialise edges and adjust p->v[v].hit.x accordingly.
  // create onb around e->omega (straight connection)
  float onb_u[3], onb_v[3];
  get_onb(p->e[v].omega, onb_u, onb_v);
  for(int k=0;k<3;k++)
    p->v[v].hit.x[k] = p->v[v-1].hit.x[k] * (1.0f-t) + p->v[v+1].hit.x[k] * t + 
      s*r*(sin_phi * onb_u[k] + cos_phi * onb_v[k]);
  // also make it a point in the volume, shader_prepare will pick up the rest below
  p->v[v].hit.prim = INVALID_PRIMID;
  p->v[v].hit.shader = p->e[v].vol.shader;

  // recreate omega/dist on the edges
  for(int k=0;k<3;k++)
    p->e[v].omega[k] = p->v[v].hit.x[k] - p->v[v-1].hit.x[k];
  p->e[v].dist = sqrtf(dotproduct(p->e[v].omega, p->e[v].omega));
  for(int k=0;k<3;k++) p->e[v].omega[k] *= 1./p->e[v].dist;
  for(int k=0;k<3;k++)
    p->e[v+1].omega[k] = p->v[v+1].hit.x[k] - p->v[v].hit.x[k];
  p->e[v+1].dist = sqrtf(dotproduct(p->e[v+1].omega, p->e[v+1].omega));
  for(int k=0;k<3;k++) p->e[v+1].omega[k] *= 1./p->e[v+1].dist;

  // test visibility of both new segments:
  if(!path_visible(p, v  )) goto fail;
  if(!path_visible(p, v+1)) goto fail;

  shader_prepare(p, v);
  shader_prepare(p, v+1);
  // sampled point so far at the rims of the volume that we fell
  // off it trying to determine the emission
  if((p->v[v].mode & s_volume) && mf_all(mf_eq(p->v[v].interior.mu_t, mf_set1(0.0f))))
    goto fail;

  // XXX hack works for lt only: reconnect camera:
  p->length = v+2; // let cam know what's the last vertex
  // we sampled this already during nee. in particular, don't use any more
  // random numbers on the last vertex (there aren't any, nee has to be
  // performed one vertex back)
  p->sensor.aperture_set = 1;
  if(p->v[0].mode & s_emit)
    edf = view_cam_connect(p);
  // TODO: also reconnect to light (in case of cosine power EDF)

  // compute transmittance and egde contribution
  const mf_t transmittance =
    shader_vol_transmittance(p, v+1) *
    shader_vol_transmittance(p, v);

#if 1 // patented pre-cancelled jacobian
  const double cosh = omega[0];
  // const double theta = acos(CLAMP(cosh, -1.0, 1.0));
  const double sinh = sqrt(MAX(0.0, 1.0-cosh*cosh));
  // compute a few needed intermediates (spherical coordinates)
  const double sp_r     = sqrt(fmax(0.0, s*s * (t-0.5)*(t-0.5) + s*s*r*r));
  const double sp_theta = acos(CLAMP((t-0.5) * s / sp_r, -1.0, 1.0));
  const double GJ       = fabs(4.0*sp_r*sin(sp_theta) / MAX(1e-7, fabs((4.0*sp_r*sp_r*cos(2.0*sp_theta) - s*s)*sinh)));

  const double f = shader_brdf(p, v-1) * transmittance * edf * p->v[v].interior.mu_s *
    path_lambert(p, v-1, p->e[v].omega) * path_lambert(p, v+1, p->e[v+1].omega);
  const double throughput = f * GJ;
  p->v[v].throughput = throughput * p->v[v-1].throughput;
#endif

  // compute jacobian as in unit test
#if 0
  const double cosh = omega[0];
  const double theta = acos(CLAMP(cosh, -1.0, 1.0));
  const double pdf_theta = hg_pdf * sin(theta) * 2.0*M_PI; // TODO: this disregards half hemisphere!

  // XXX TODO: can't compute this case, at least not without cancelling!
  if(theta < 1e-6) goto fail;

  // const double sinh2 = 1.0 - cosh*cosh;
  // const double sinh = sqrt(sinh2);
  // sample t and compute r
  // const double s = 1.0; // global scale factor. s=1 means radius of sphere = 1/2 and distance x0--x2=1.
  // const double t = t_sample(cosh, r2);
  // const double r = sqrt(1.0/(4.0*sinh2) - (0.5-t)*(0.5-t)) - sqrt(1.0/(4.0*sinh2) - 0.25);
  // const double phi = 2.0 * M_PI * r1;
  // const double x1[3] = {s*t, s*r * sin(phi), s*r * cos(phi) };


  // compute a few needed intermediates (spherical coordinates)
  const double sp_r     = sqrt(fmax(0.0, s*s * (t-0.5)*(t-0.5) + s*s*r*r));
  const double sp_theta = acos(CLAMP((t-0.5) * s / sp_r, -1.0, 1.0));
  //             sp_phi = phi; // <= this one stays, yay.
  // assert(theta == theta);
  // assert(sp_theta == sp_theta);
  // assert(sp_r == sp_r);
  // test_float(t, (sp_r / s * cos(sp_theta) + 0.5));
  // // both equations are correct:
  // test_float(r2, (1.0/(2.0*theta)*(theta + asin((2.0*t-1.0)*sin(theta)))));
  // test_float(r2, (1.0/(2.0*theta)*(theta - asin((1.0-2.0*t)*sin(theta)))));

  // // test side lengths/law of cosines
  // const double d_1_sp = sqrt(sp_r*sp_r + s*s/4.0 + sp_r*s*cos(sp_theta));
  // const double d_1_cy = sqrt(s*s*t*t + s*s*r*r);
  // // test_float(d_1_sp, d_1_cy);

  // const double d_2_sp = sqrt(sp_r*sp_r + s*s/4.0 - sp_r*s*cos(sp_theta));
  // const double d_2_cy = sqrt(s*s*(1.0-t)*(1.0-t) + s*s*r*r);
  // // test_float(d_2_sp, d_2_cy);

#if 0
  // compute the vertex both ways and make sure it's close enough
  const double x1c[3] = {sp_r * cos(sp_theta) + 0.5*s, sp_r*sin(sp_theta) * sin(phi), sp_r*sin(sp_theta) * cos(phi)};
  test_float(x1c[0], x1[0]);
  test_float(x1c[1], x1[1]);
  test_float(x1c[2], x1[2]);
#endif
  // const double delta = 1e-7;

  // TODO: cancel this "den" because it blows up for theta near 0
  // compute Jacobian from t,theta to spherical coordinates sp_r,sp_theta
  const double den = fmax(1e-6, 16.0*pow(sp_r, 4.0) + pow(s, 4.0) - 8.0 *sp_r*sp_r *s*s * cos(2*sp_theta));
  const double DthetaDsp_r     = 4.0*s*(4.0*sp_r*sp_r + s*s)*sin(sp_theta) / den;
  const double DthetaDsp_theta = 4.0*sp_r*s*(s*s - 4.0*sp_r*sp_r)*cos(sp_theta) / den;

  assert(den==den);
  assert(DthetaDsp_theta == DthetaDsp_theta);
  assert(DthetaDsp_r     == DthetaDsp_r);
#if 0
  // const double test_theta2 = M_PI - acos((d_1_cy*d_1_cy + d_2_cy*d_2_cy - s*s)/(2.0*d_1_cy*d_2_cy));
  // const double test_theta3 = M_PI - acos((sp_r*sp_r - s*s/4.0)/(d_1_cy*d_2_cy));
  const double test_theta = theta_from_sp(s, sp_r, sp_theta);
  test_float(test_theta, theta);
  const double DthetaDsp_r_fd     = (theta_from_sp(s, sp_r+delta, sp_theta) - theta_from_sp(s, sp_r-delta, sp_theta))/(2.0*delta);
  const double DthetaDsp_theta_fd = (theta_from_sp(s, sp_r, sp_theta+delta) - theta_from_sp(s, sp_r, sp_theta-delta))/(2.0*delta);
  test_float(DthetaDsp_r,     DthetaDsp_r_fd);
  test_float(DthetaDsp_theta, DthetaDsp_theta_fd);
#endif

  // outer derivative:
  const double DxiDt     = sin(theta) / (theta * sqrt(fmax(0.0, 1.0 - (2.0*t-1.0)*(2.0*t-1.0)*sin(theta)*sin(theta))));
  const double DxiDtheta = 1.0/(2.0*theta*theta) * (asin(CLAMP((1.0-2.0*t)*sin(theta), -1.0, 1.0)) +
                           (2.0*t-1.0)*theta*cos(theta)/sqrt(fmax(0.0, 1.0 - (1.0-2.0*t)*(1.0-2.0*t)*sin(theta)*sin(theta))));
  assert(DxiDt==DxiDt);
  assert(DxiDtheta==DxiDtheta);
#if 0 // seems to work, good:
  const double test_xi = xi_from_ttheta(s, t, theta);
  test_float(test_xi, r2);
  const double DxiDt_fd     = (xi_from_ttheta(s, t+delta, theta) - xi_from_ttheta(s, t-delta, theta))/(2.0*delta);
  const double DxiDtheta_fd = (xi_from_ttheta(s, t, theta+delta) - xi_from_ttheta(s, t, theta-delta))/(2.0*delta);
  test_float(DxiDt,      DxiDt_fd    );
  test_float(DxiDtheta,  DxiDtheta_fd);
#endif

  // additional inner derivatives:
  const double DtDsp_theta = - sp_r/s * sin(sp_theta);
  const double DtDsp_r     =    1.0/s * cos(sp_theta);

  assert(DtDsp_theta==DtDsp_theta);
  assert(DtDsp_r==DtDsp_r);
#if 0
  const double test_t = t_from_sp(s, sp_r, sp_theta);
  test_float(test_t, t);
  const double DtDsp_r_fd     = (t_from_sp(s, sp_r+delta, sp_theta) - t_from_sp(s, sp_r-delta, sp_theta))/(2.0*delta);
  const double DtDsp_theta_fd = (t_from_sp(s, sp_r, sp_theta+delta) - t_from_sp(s, sp_r, sp_theta-delta))/(2.0*delta);
  test_float(DtDsp_r,     DtDsp_r_fd    );
  test_float(DtDsp_theta, DtDsp_theta_fd);
#endif

  // finally the 2D chain rule:
  const double DxiDsp_r     = DxiDtheta * DthetaDsp_r     + DxiDt * DtDsp_r;
  const double DxiDsp_theta = DxiDtheta * DthetaDsp_theta + DxiDt * DtDsp_theta;

  assert(DxiDsp_theta==DxiDsp_theta);
  assert(DxiDsp_r==DxiDsp_r);

#if 0 // test derivatives against finite differences:
  const double xi = xi_from_sp(s, sp_r, sp_theta);
  test_float(xi, r2);
  const double DxiDsp_r_fd     = (xi_from_sp(s, sp_r + delta, sp_theta) - xi_from_sp(s, sp_r - delta, sp_theta))/(2.0*delta);
  const double DxiDsp_theta_fd = (xi_from_sp(s, sp_r, sp_theta + delta) - xi_from_sp(s, sp_r, sp_theta - delta))/(2.0*delta);
  test_float(DxiDsp_r,     DxiDsp_r_fd);
  test_float(DxiDsp_theta, DxiDsp_theta_fd);
#endif

  // now the jacobian itself.
  const double J = 1.0/fabs(DthetaDsp_r * DxiDsp_theta - DthetaDsp_theta * DxiDsp_r);
  // const double J = 1.0/fabs(DthetaDsp_r_fd * DxiDsp_theta_fd - DthetaDsp_theta_fd * DxiDsp_r_fd);

  assert(J==J);

  // pdf(r2)    = 1
  // pdf(theta) = pdf_phase_func_solid_angle * sinh *2*pi  (because we are not sampling phi with this)
  // pdf(phi)   = 1/2pi
  const double pdf_phi = 1.0/(2.0*M_PI);
  const double our_pdf = pdf_phi * pdf_theta;

  const double f = shader_brdf(p, v-1) * transmittance * edf * p->v[v].interior.mu_s * path_G(p, v) * path_G(p, v+1);
  const double throughput = f / our_pdf * J * sin(sp_theta) * sp_r * sp_r;
  p->v[v].throughput = throughput * p->v[v-1].throughput;
#endif

  // XXX can we explicitly cancel the g terms?
#if 0
  const float theta = acosf(CLAMP(omega[0],-1,1));
  // c = G * G / pdf
  const float sinh = sqrtf(sinh2);
  const float J = theta / (2.0f * sinh2 * sinh) * (cosf((2.0*dist_r-1.0)*theta) - cosf(theta));
  const float pdf_x1 = /* cancelled: hg_pdf */  /* p(xi==dist_r) = 1 */ sinh / J;
  const float c = s* // vertex area pdf is actually lower if we scale up the volume
    2.0*M_PI* // sampling phi uniformly on the circle
    // path_lambert(p, v-1, p->e[v].omega) * path_lambert(p, v+1, p->e[v+1].omega)
    path_G(p, v) * path_G(p, v+1)
    / pdf_x1;
  // recompute bsdf based on new direction
  bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1
  // T*T*f_r * EDF * mu_s * c
  p->v[v].throughput = mf_mul(mf_mul(mf_mul(p->v[v-1].throughput, bsdf), mf_mul(transmittance, edf)), mf_mul(p->v[v].interior.mu_s, mf_set1(c)));
#endif

#if 0
  // compute and multiply:
  // cos0 cos2
  // mu_s(v)    // potentially spectral
  // 2pi s^3    // i think all of these cancel ( p(phi) = 2pi and we blow up normalised vertex area measure by scaling it back)
  // theta/sin(theta)
  const float theta = acosf(CLAMP(omega[0],-1,1));
  // fudge factor in front makes it a tad too dark.
  // at about double distance this error halves :(
  // d = 175 => 24%
  // d = 345 => 12%
  // totally *not* another dependency on s though.
  float c = // 1.0/(M_PI*s) *
    1.0/(s*s) *
    theta/sinf(theta) *
    // theta *
    path_lambert(p, v-1, p->e[v].omega) * path_lambert(p, v+1, p->e[v+1].omega);
  // recompute bsdf based on new direction
  bsdf = shader_brdf(p, v-1); // also set mode on vertex v-1
  p->v[v].throughput = mf_mul(mf_mul(mf_mul(p->v[v-1].throughput, bsdf), mf_mul(transmittance, edf)), mf_mul(p->v[v].interior.mu_s, mf_set1(c)));
#endif
#if 0
  // XXX these two do not match!!
  // XXX the above is off by a factor of two and this is just nuts:
  // mean image brightness here seems okayish but has apparently a lot of variance.
  // this is p(t) * p(theta) * p(phi) * sin(theta) (last term is double in t_pdf and hg_pdf, so we cancel it out)
  // also we have to scale the pdf to match the normalisation of the scattering geometry to [0,1], so s^3
  // or no s^3? maybe it cancels mysteriously with something.
  const double pdf_t_theta_phi = t_pdf(t, omega[0]);// / (sin(theta));// * (double)s * (double)s * (double)s);
  const double eval_bsdf = shader_brdf(p, v-1) *
    shader_brdf(p, v); // this is bsdf * mu_s * phase func
  // this is the throughput to be multiplied to p->v[v-1].throughput (whatever was incoming to the single scattering point)
  double throughput =
    // 1.0/(s*s)*

    s*(double)s*s*(double)s *
    1.0/((double)p->e[v].dist*(double)p->e[v].dist*(double)p->e[v+1].dist*(double)p->e[v+1].dist) * // dist term of the two G
    // XXX CLAMP(1.0/((double)p->e[v].dist*(double)p->e[v].dist*(double)p->e[v+1].dist*(double)p->e[v+1].dist), 0, 100) * // dist term of the two G
    // path_lambert(p, v-1, p->e[v].omega) * path_lambert(p, v+1, p->e[v+1].omega) *
    // fabsf(dotproduct(p->v[v-1].hit.n, p->e[v].omega)) *       // cosine at single scattering
    // fabsf(dotproduct(p->v[v+1].hit.n, p->e[v+1].omega)) *     // cosine at camera
    // eval_bsdf *
    // transmittance * edf *
    // 1.0 / hg_pdf *
    1.0 / pdf_t_theta_phi;
  // these look good, we're mostly sampling towards the ends of the NEE segment
  // fprintf(stderr, "dists s %g d1 %g d2 %g cos theta %g\n", s, p->e[v].dist, p->e[v+1].dist, omega[0]);
  // fprintf(stderr, "terms are pdf %g bsdf %g G %g cos %g T %g edf %g res %g\n",
  //     pdf_t_theta_phi, eval_bsdf, 
  //   1.0/((double)p->e[v].dist*(double)p->e[v].dist*(double)p->e[v+1].dist*(double)p->e[v+1].dist),
  //   fabsf(dotproduct(p->v[v-1].hit.n, p->e[v].omega)) *
  //   fabsf(dotproduct(p->v[v+1].hit.n, p->e[v+1].omega)),
  //   transmittance, edf, throughput);

  // throughput = shader_brdf(p, v)/ hg_pdf;
  fprintf(stderr, "throughput %g\n", throughput);
  // mf_t throughput = mf_div(mf_mul(
  //     mf_mul(p->v[v-1].throughput, mf_set1((float)(path_G(p, v) * path_G(p, v+1)))),
  //     mf_mul(mf_mul(bsdf, transmittance), mf_mul(transmittance, edf))),
  //     pdf);
  // double thr_opt = bsdf * transmittance * edf * p->v[v].interior.mu_s * c;
  // const double rtio = thr_opt / throughput;
  // fprintf(stderr, "thr E %g f/p %g ratio %g s %g\n", thr_opt, throughput, rtio, s);
  p->v[v].throughput = p->v[v-1].throughput * throughput; // XXX see what the complicated eval says
#endif

  p->v[v+1].throughput = p->v[v].throughput; // just set both to the same thing

  p->throughput = p->v[v+1].throughput; // already contains light/sensor edf

  if(0)
  {
fail:
    p->v[v].pdf = mf_set1(0.0f);
    p->throughput = mf_set1(0.0f);
    p->v[v].throughput = mf_set1(0.0f);
    p->v[v].flags = s_none;
    p->v[v].mode = s_absorb;
    p->v[v].rand_cnt = s_dim_num_nee;
    p->v[v+1] = p->v[v];
    // need two more, because pop will also pop two
    p->length = v+2; // constructed vertex v (even if it absorbs)
    return 0;
  }
  p->v[v  ].rand_cnt = s_dim_num_nee;
  p->v[v+1].rand_cnt = s_dim_num_nee;
  p->length = v+2; // constructed vertex v and v+1

  const mf_t pdf_nee = p->v[v+1].pdf;
  p->v[v+1].total_throughput = p->v[v+1].throughput; // store for path_pop()
  p->v[v+1].pdf = pdf_nee;
  // XXX all wrong! (but currently unused)
  p->v[v].pdf = mf_set1(1.0f/s*t_pdf(t, omega[0]) * hg_pdf);

  // TODO: validate our estimator by evaluating throughput as full measurement / this:
#if 0 // DEBUG
  // parts of the measurement we care about: G1 G2 fs
  // float fs = sample_eval_hg_fwd(p->v[v].interior.mean_cos, p->e[v].omega, p->e[v+1].omega);   // should probably eval the phase function via shader_brdf()
  // float lam = path_lambert(p, v-1, p->e[v].omega) * path_lambert(p, v+1, p->e[v+1].omega);
  // float measurement = path_G(p, v) * path_G(p, v+1) / lam;// * fs;
  double measurement = 1./(p->e[v].dist*p->e[v].dist) *
                       1./(p->e[v+1].dist*p->e[v+1].dist);
  // fprintf(stderr, "dots %g %g\n", omega[0], dotproduct(p->e[v].omega, p->e[v+1].omega));
  // float mvpdf = 2*M_PI/s*t_pdf(t, omega[0]);// * hg_pdf;
  // XXX there is an s^2 here that we can't have in the estimator, mean img brightness will be way wrong!
  double mvpdf = theta * 2.0*M_PI/(s*s*s)*t_pdf(t, omega[0]);// * hg_pdf;
  fprintf(stderr, "ratio %g\n", measurement / mvpdf);
#endif

  return 0;
}
