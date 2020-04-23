/*
    This file is part of corona-13.
    copyright (c) 2004--2017 johannes hanika

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
#include "camera.h"
#include "view.h"
#include "sampler_common.h"
#include "pathspace.h"
#include "pointsampler.h"
#include "spectrum.h"
#include "display.h"
#include "aperture.h"
#include "lens/raytrace.h"
#include "lens/sampling.h"

#include <math.h>
#include <unistd.h>
#define CAMERA_BLADES 9

static lens_t lens;

int camera_global_init()
{
  const char *lensname = CAMERA_LENS;
  char exe[512];
  readlink("/proc/self/exe", exe, sizeof(exe));
  char *c = exe + strlen(exe);
  for(;*c != '/' && c>exe;c--);
  snprintf(c, sizeof(exe) - (c - exe), "/camera/%s/table", lensname);
  int err = lens_read(&lens, exe);
  if(err)
  {
    fprintf(stderr, "[lens] failed to read `%s'!\n", exe);
    exit(1);
  }
  fprintf(stdout, "[lens] loading `%s'\n", lensname);
  lens_canonicalize_name(lensname, lens.name);
  return err;
}

void camera_display_info(const camera_t *c)
{
  if(c->exposure_value > 6)
    display_print(rt.display, 0, 0, " %s   1/%.0f f/%.1f iso%d", lens.name, roundf(1.0/view_tv2time(c->exposure_value)),
        view_av2fstop(c->aperture_value), (int)c->iso);
  else
    display_print(rt.display, 0, 0, " %s   %.1f\" f/%.1f iso%d", lens.name, view_tv2time(c->exposure_value),
        view_av2fstop(c->aperture_value), (int)c->iso);
}

void camera_print_info(const camera_t *c, FILE *f)
{
  fprintf(f, "camera   : ray traced optics lens `%s'\n", lens.name);
  fprintf(f, "  focus  : %fm\n", c->focus/10.0f);
  fprintf(f, "  sensor : %fmm (offset)\n", c->focus_sensor_offset);
  fprintf(f, "  film   : %dmm x %dmm\n", (int)(c->film_width*100.0), (int)(c->film_height*100.0 + 0.5f));
  fprintf(f, "  iso    : %d\n", (int)c->iso);
  const float time = view_tv2time(c->exposure_value);
  if(time < 0.5)
    fprintf(f, "  Tv     : 1/%.0f\n", roundf(1.0/time));
  else
    fprintf(f, "  Tv     : %.1f''\n", time);
  fprintf(f, "  Av     : f/%.1f\n", view_av2fstop(c->aperture_value));
}

static inline float camera_aperture_radius(const camera_t *c)
{
  return MIN(lens.aperture_housing_radius, lens.focal_length/(2.0f * view_av2fstop(c->aperture_value)));
}

void camera_set_focus(camera_t *c, float dist)
{
  c->focus = dist;
  const float dm2mm = 100.0f;
  // camera space vector to v1:
  const float target[3] = { 0.0, 0.0, dm2mm * dist };
  // trace a couple of adjoint rays from there to the sensor and see where we need to put the sensor plane.
  float ap[3] = {0.0f}; // aperture position
  float spos[3], sdir[3];
  float opos[3], odir[3];
  float off = 0.0f;
  int cnt = 0;
  const float arad = camera_aperture_radius(c);
  const int S = 4;
  for(int s=1;s<=S;s++)
  for(int k=0;k<2;k++)
  {
    ap[0] = ap[1] = 0.0f;
    ap[k] = arad * s/(S+1.0f);
    float transmittance = lens_lt_sample_aperture(
        &lens, 550.0f, target, ap, spos, sdir, opos, odir, 0);
    // fprintf(stderr, "%d sensor tr %g pos/dir %g %g %g -- %g %g %g\n",
    //     k, transmittance,
    //     spos[0], spos[1], spos[2],
    //     sdir[0], sdir[1], sdir[2]);
    if(transmittance > 0 && -sdir[k] > 0)
    {
      off += -spos[k]/sdir[k];
      cnt ++;
      // fprintf(stderr, "estimate %g %d => %g\n", off, cnt, off/cnt);
    }
  }
  off /= cnt; // average guesses
  // the focus plane/sensor offset:
  fprintf(stderr, "[camera] focus shift  : %.4fmm (sensor offset)\n", off);
  // negative because of reverse direction
  if(off == off)
  {
    const float limit = 45.0f;
    if(off > limit) c->focus_sensor_offset = limit;
    else if(off < -limit) c->focus_sensor_offset = -limit;
    else c->focus_sensor_offset = off; // in mm
  }
}

float camera_sample(const camera_t *c, path_t *p)
{
  const float i = p->sensor.pixel_set ? (p->sensor.pixel_i / view_width())  : pointsampler(p, s_dim_image_x);
  const float j = p->sensor.pixel_set ? (p->sensor.pixel_j / view_height()) : pointsampler(p, s_dim_image_y);
  const float r1 = p->sensor.aperture_set ? p->sensor.aperture_x : pointsampler(p, s_dim_aperture_x);
  const float r2 = p->sensor.aperture_set ? p->sensor.aperture_y : pointsampler(p, s_dim_aperture_y);

  const float dm2mm = 100.0f;

  float spos[3] = {0}, sdir[3] = {0}, apos[3] = {0},
        opos[3] = {0}, odir[3] = {0};
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);
  lens_aperture_sample(apos+0, apos+1, r1, r2, aperture_radius, num_blades);
  const float A = lens_aperture_area(aperture_radius, num_blades);

  // sample point on the sensor in mm:
  spos[0] = dm2mm*(i-.5f)*c->film_width;
  spos[1] = dm2mm*(j-.5f)*c->film_height;
  spos[2] = -c->focus_sensor_offset;

  float n[3] = {0};

  // get direction via newton iteration:
  float transmittance = lens_pt_sample_aperture(
      &lens, p->lambda,
      spos, sdir, apos, opos, odir, n);
  if(transmittance <= 0.0f) return 0.0f;

  // remember pixel coordinates:
  p->sensor.pixel_i = (spos[0]/(dm2mm*c->film_width)  + .5f)*view_width();
  p->sensor.pixel_j = (spos[1]/(dm2mm*c->film_height) + .5f)*view_height();
  p->sensor.pixel_dir_i = sdir[0]/sdir[2];
  p->sensor.pixel_dir_j = sdir[1]/sdir[2];
  p->sensor.aperture_x = r1;
  p->sensor.aperture_y = r2;

  for(int k=0;k<3;k++) opos[k] /= dm2mm;

  // fprintf(stderr, "opos %g %g %g\n", opos[0], opos[1], opos[2]);
  view_cam_init_frame(p, &p->v[0].hit);
  for(int k=0;k<3;k++)
  {
    p->v[0].hit.x[k] += opos[0] * p->v[0].hit.a[k] + opos[1] * p->v[0].hit.b[k] + opos[2] * p->v[0].hit.n[k];
    p->e[1].omega[k] = odir[0] * p->v[0].hit.a[k] + odir[1] * p->v[0].hit.b[k] + odir[2] * p->v[0].hit.n[k];
  }

  for(int k=0;k<3;k++) p->v[0].hit.n[k] = n[k];
  // TODO: recompute frame with a/b, too?

  // also clip to valid pixel range.
  if(p->sensor.pixel_i < 0.0f || p->sensor.pixel_i >= view_width() ||
     p->sensor.pixel_j < 0.0f || p->sensor.pixel_j >= view_height())
  {
    p->v[0].material_modes = p->v[0].mode = s_absorb;
    return 0.0f;
  }

  // compute geometric term and pdfs
  const float inv_p = c->film_width*c->film_height*A/(dm2mm*dm2mm); // XXX blows up?
  const float sensor = 100.0f*view_tv2time(c->exposure_value);

  // adjust pathspace parameters on vertices:
  p->v[0].hit.prim = INVALID_PRIMID;
  p->v[0].hit.shader = -1;
  p->v[0].pdf = 1.0f; // first two vertices created in one go, pdf is only in v[1]
  p->v[0].material_modes = p->v[0].mode = s_sensor;
  p->v[0].rand_cnt = s_dim_num_pt_beg;

  p->v[1].rand_beg = p->v[0].rand_beg + p->v[0].rand_cnt;
  p->v[1].rand_cnt = 1; // only free path length
  p->v[1].pdf = 1.0f/inv_p; // will be transformed by geometric term later on.

#if 0
  // XXX compute this as G*|Tp| and pass out of something?
  // cos^2 theta(sen) twice (once for propagation along optical axis and once for angles)
  const float det = 1.0f/(1.0f + sen[2]*sen[2] + sen[3]*sen[3]);
  // the determinant sensor to aperture is required both for sample and connect, since
  // both now sample a point on the aperture.
  // this transforms position on aperture to projected solid angle direction on sensor:
  // path tracing needs it to convert pdfs
  const float deta = lens_det_aperture_to_sensor(sen, c->focus_sensor_offset) * det*det;
#endif
  return transmittance * sensor;//*inv_p;
}

float camera_pdf(const camera_t *c, const path_t *p, int v)
{
  const float dm2mm = 100.0f;
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);
  // TODO: sure this is correct? do we need vertex area measure of point on outgoing pupil? (same as camera_sample..)
  // TODO: anyways should be made consistent with _connect which seems to work on outer pupil vertex area measure.
  const float A = lens_aperture_area(aperture_radius, num_blades)/(dm2mm*dm2mm);
  return 1.0f/(c->film_width*c->film_height*A);
}

// TODO: port to manifold thing
float camera_connect(const camera_t *c, path_t *p)
{
  return 0.0f; // XXX
#if 0
  //             path tracing case               mutation case, end vertex there already              light tracing case, create a new vertex at the end     
  const int v0 = (p->v[0].mode & s_sensor) ? 0 : ((p->v[p->length-1].mode & s_sensor) ? p->length-1 : p->length);
  const int v1 = (p->v[0].mode & s_sensor) ? 1 : ((p->v[p->length-1].mode & s_sensor) ? p->length-2 : p->length-1);
  const int e  = (p->v[0].mode & s_sensor) ? 1 : ((p->v[p->length-1].mode & s_sensor) ? p->length-1 : p->length);

  const float ox = p->sensor.aperture_set ? p->sensor.aperture_x : pointsampler(p, s_dim_nee_x);
  const float oy = p->sensor.aperture_set ? p->sensor.aperture_y : pointsampler(p, s_dim_nee_y);
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);

  const float dm2mm = 100.0f;

  // [x,y,dx,dy,lambda] on the sensor, aperture, and outgoing pupil:
  float sen[5], ape[5], out[5];
  sen[4] = ape[4] = p->lambda/1000.0f; // nanometers to micrometers

  view_cam_init_frame(p, &p->v[v0].hit);
  const float R = lens_outer_pupil_curvature_radius;
#ifdef LENS_OUTER_PUPIL_LT
  out[0] = cosf(2*M_PI*ox)*sqrtf(oy)*lens_outer_pupil_radius;
  out[1] = sinf(2*M_PI*ox)*sqrtf(oy)*lens_outer_pupil_radius;
  const float out2 = sqrtf(R*R - out[0]*out[0] - out[1]*out[1]);
  const float cos_theta = out2/fabsf(R); // p(out) = 1/cos_theta
  // if(cos_theta < lens_field_of_view) return 0.0f;

  float campos[3] = {out[0], out[1], out2-R};
  float view[3];
  for(int k=0;k<3;k++)
  {
    p->v[v0].hit.x[k] += p->v[v0].hit.a[k] * campos[0] / dm2mm + p->v[v0].hit.b[k] * campos[1] / dm2mm + p->v[v0].hit.n[k] * campos[2] / dm2mm;
    view[k] = p->v[v1].hit.x[k] - p->v[v0].hit.x[k];
  }
  p->e[e].dist = sqrtf(dotproduct(view, view));
  if(p->v[0].mode & s_sensor)
    for(int k=0;k<3;k++) p->e[e].omega[k] = view[k] / p->e[e].dist;
  else
    for(int k=0;k<3;k++) p->e[e].omega[k] = - view[k] / p->e[e].dist;
  float camdir[3] = {dotproduct(p->e[e].omega, p->v[v0].hit.a), dotproduct(p->e[e].omega, p->v[v0].hit.b), dotproduct(p->e[e].omega, p->v[v0].hit.n)};
  // XXX recompute full frame:
  float n[3] = {0.0f};
  for(int k=0;k<3;k++)
    n[k] += p->v[v0].hit.a[k] * out[0]/R + p->v[v0].hit.b[k] * out[1]/R + p->v[v0].hit.n[k] * cos_theta;
  for(int k=0;k<3;k++) p->v[v0].hit.n[k] = n[k];
  if(p->v[0].mode & s_emit) for(int k=0;k<3;k++) camdir[k] = - camdir[k];
  lens_csToSphere(campos, camdir, out, out+2, -R, R);

  // now iterate to find this outgoing configuration:
  sen[0] = sen[1] = sen[2] = sen[3] = 0.0f;
  float outp[5] = {0.0f}, J[25] = {0.0f}, Ji[25] = {0.0f}, outdiff[5] = {0.0f};
  float transmittance = 0.0f;
  float olderr = FLT_MAX;
  for(int it=0;it<100;it++)
  {
    transmittance = lens_evaluate(sen, outp);
    for(int k=0;k<4;k++) outdiff[k] = out[k] - outp[k];
    const float err = outdiff[0]*outdiff[0] + outdiff[1]*outdiff[1] + outdiff[2]*outdiff[2] + outdiff[3]*outdiff[3];
    if(err < 1e-4f) break;
    transmittance = 0.0f;
    if(it >= 10 && err > olderr) break;
    lens_evaluate_jacobian(sen, J);
    lens_mat5inv( (const float (*)[5])J, (float (*)[5])Ji);
    float sendiff[5] = {0.0f};
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
        sendiff[j] += Ji[5*j+i]*outdiff[i];

    for(int j=0;j<4;j++) sen[j] += sendiff[j];
    olderr = err;
  }
  if(transmittance <= 0.0f) return 0.0f;

  lens_evaluate_aperture(sen, ape);
  if(!lens_clip_aperture(ape[0], ape[1], aperture_radius, num_blades)) return 0.0f;

  // density is on spherical outer pupil surface (cos/disk area)
  // this does not match camera_connect(), since this is just for debugging purposes (don't use it for real!)
  p->v[v0].pdf = dm2mm*dm2mm*cos_theta/(M_PI*lens_outer_pupil_radius*lens_outer_pupil_radius);
#else // sample point on aperture:
  lens_sample_aperture(ape+0, ape+1, ox, oy, aperture_radius, num_blades);

  float totarget[3];
  for(int k=0;k<3;k++) totarget[k] = p->v[v1].hit.x[k] - p->v[v0].hit.x[k];
  // camera space vector to v1:
  const float target[3] =
  {
    dm2mm * dotproduct(p->v[v0].hit.a, totarget),
    dm2mm * dotproduct(p->v[v0].hit.b, totarget),
    dm2mm * dotproduct(p->v[v0].hit.n, totarget)
  };

  float transmittance = lens_lt_sample_aperture(target, ape, sen, out, p->lambda/1000.0f);
  if(transmittance <= 0.0f) return 0.0f;
  float campos[3], camdir[3];
  lens_sphereToCs(out, out+2, campos, camdir, -R, R);
  for(int k=0;k<3;k++)
  {
    p->v[v0].hit.x[k] += p->v[v0].hit.a[k] * campos[0] / dm2mm + p->v[v0].hit.b[k] * campos[1] / dm2mm + p->v[v0].hit.n[k] * campos[2] / dm2mm;
    p->e[e].omega[k] = p->v[v1].hit.x[k] - p->v[v0].hit.x[k];
  }
  p->e[e].dist = sqrtf(dotproduct(p->e[e].omega, p->e[e].omega));
  if(p->v[0].mode & s_sensor)
    for(int k=0;k<3;k++) p->e[e].omega[k] /= p->e[e].dist;
  else
    for(int k=0;k<3;k++) p->e[e].omega[k] /= - p->e[e].dist;
  // XXX recompute full frame:
  float n[3] = {0.0f};
  for(int k=0;k<3;k++)
    n[k] += p->v[v0].hit.a[k] * campos[0]/R + p->v[v0].hit.b[k] * campos[1]/R + p->v[v0].hit.n[k] * (campos[2]+R)/fabsf(R);
  for(int k=0;k<3;k++) p->v[v0].hit.n[k] = n[k];

  // this transforms position on aperture to plane/plane direction on sensor.
  // light tracing needs it as an intermediate step to convert sampling density from point on aperture to point on outer pupil.
  const float deta = lens_det_aperture_to_sensor(sen, 0);
  // converts sensor plane/plane direction to vertex area measure on outer pupil.
  // we transform pdf from vertex area on aperture -> direction on sensor -> position on outer pupil
  const float deto = lens_det_sensor_to_outer_pupil(sen, out, 0);
  const float det_ao = deta * deto;
#if 1 // be consistent with mis:
  p->v[v0].pdf = dm2mm*dm2mm/(lens_aperture_area(aperture_radius, num_blades) * det_ao);
#else
  transmittance *= det_ao;
  p->v[v0].pdf = dm2mm*dm2mm/lens_aperture_area(aperture_radius, num_blades);
#endif

  // XXX DEBUG mul to pt instead?
  // transmittance /= fabsf(dotproduct(p->e[e].omega, p->v[v0].hit.n));
#endif
  // sensor pos - sensor direction because the directions are always going to the outside
  p->sensor.pixel_i = .5f*view_width()  + (sen[0]-sen[2]*c->focus_sensor_offset)/(dm2mm*c->film_width) * view_width();
  p->sensor.pixel_j = .5f*view_height() + (sen[1]-sen[3]*c->focus_sensor_offset)/(dm2mm*c->film_height)* view_height();
  p->sensor.pixel_dir_i = sen[2];
  p->sensor.pixel_dir_j = sen[3];
  p->sensor.aperture_x = ox;
  p->sensor.aperture_y = oy;
  const float sensor = 100.0f*view_tv2time(c->exposure_value);
  p->v[v0].material_modes =
  p->v[v0].mode = s_sensor;
  p->v[v0].flags = s_none;
  p->v[v0].hit.prim = INVALID_PRIMID;
  p->v[v0].hit.shader = -1;
  p->v[v0].shading.roughness = 1.0f; // instruct path space that we don't need an additional cosine

  // crop at inward facing pupil:
  const float px = sen[0] + sen[2] * lens_focal_length, py = sen[1] + sen[3]*lens_focal_length;
  if(px*px + py*py > lens_inner_pupil_radius*lens_inner_pupil_radius) return 0.0f;
  if(!(p->sensor.pixel_i >= 0 && p->sensor.pixel_j >= 0 && p->sensor.pixel_i < view_width() && p->sensor.pixel_j < view_height())) return 0.0f;
  return transmittance * sensor/p->v[v0].pdf;
#endif
}

float camera_pdf_connect(const camera_t *c, const path_t *p, int v)
{
  return 0.0f;// XXX
#if 0
#if 1 // needs to be in sync with above
  const int v1 = v == 0 ? 1 : v-1;

  const float ox = p->sensor.aperture_x;
  const float oy = p->sensor.aperture_y;
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);

  const float dm2mm = 100.0f;

  // [x,y,dx,dy,lambda] on the sensor, aperture, and outgoing pupil:
  float sen[5], ape[5], out[5];
  sen[4] = ape[4] = p->lambda/1000.0f; // nanometers to micrometers

  hit_t frame;
  view_cam_init_frame(p, &frame);
  lens_sample_aperture(ape+0, ape+1, ox, oy, aperture_radius, num_blades);

  float totarget[3];
  for(int k=0;k<3;k++) totarget[k] = p->v[v1].hit.x[k] - frame.x[k];
  // camera space vector to v1:
  const float target[3] =
  {
    dm2mm * dotproduct(frame.a, totarget),
    dm2mm * dotproduct(frame.b, totarget),
    dm2mm * dotproduct(frame.n, totarget)
  };

  float transmittance = lens_lt_sample_aperture(target, ape, sen, out, p->lambda/1000.0f);
  if(transmittance <= 0.0f) return 0.0f;

  // this transforms position on aperture to plane/plane direction on sensor.
  // light tracing needs it as an intermediate step to convert sampling density from point on aperture to point on outer pupil.
  const float deta = lens_det_aperture_to_sensor(sen, 0);
  // converts sensor plane/plane direction to vertex area measure on outer pupil.
  // we transform pdf from vertex area on aperture -> direction on sensor -> position on outer pupil
  const float deto = lens_det_sensor_to_outer_pupil(sen, out, 0);
  const float det_ao = deta * deto;
  return dm2mm*dm2mm/(lens_aperture_area(aperture_radius, num_blades) * det_ao);
#else
  // sampled the outward facing pupil:
  const float dm2mm = 100.0f;
  const float aperture_radius = camera_aperture_radius(c);
  return dm2mm*dm2mm/lens_aperture_area(aperture_radius, CAMERA_BLADES);
#endif
#endif
}

float camera_mutate_aperture(const camera_t *c, path_t *p, const float r1, const float r2, const float step)
{
  // mutate point on aperture via kelemen mutation, call camera_sample again.
  p->sensor.aperture_x = sample_mutate_rand(p->sensor.aperture_x, r1, step);
  p->sensor.aperture_y = sample_mutate_rand(p->sensor.aperture_y, r2, step);
  p->sensor.aperture_set = 1;
  p->sensor.pixel_set = 1;
  return camera_sample(c, p);
}

float camera_pdf_mutate_aperture(const camera_t *c, const path_t *curr, const path_t *tent, const float step)
{
  // this is, again, vertex area of aperture not outgoing pupil.
  const float dm2mm = 100.0f;
  const float aperture_radius = camera_aperture_radius(c);
  return dm2mm*dm2mm/(lens_aperture_area(aperture_radius, CAMERA_BLADES) * step*step);
}

float camera_eval(const camera_t *c, const path_t *p)
{
  if(!(p->sensor.pixel_i >= 0 && p->sensor.pixel_j >= 0 && p->sensor.pixel_i < view_width() && p->sensor.pixel_j < view_height())) return 0.0f;
  const float sensor = 100.0f*view_tv2time(c->exposure_value);
  return sensor;
}
