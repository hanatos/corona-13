/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CAMERA_LENS_RT_H
#define CAMERA_LENS_RT_H

#include "camera_common.h"
#include "pathspace.h"
#include "pointsampler.h"
#include "spectrum.h"
#include "display.h"

#include "camera/init.h"

// DEBUG XXX this is just for the paper and totally useless for production use!


// fake function in micrometers:
static inline float spectrum_eta_from_abbe_um(const float n_d, const float V_d, const float lambda)
{
    float A, B;
      spectrum_cauchy_from_abbe(n_d, V_d, &A, &B);

        return A + B/(lambda*lambda);
}

// this is ugly, in the hope it won't ever compile after submission :)
#include "../optics/test/gui/lenssystem.hh" // have to make sure we link against the corresponding c++ .o as well.
#include "../optics/test/gui/raytrace.h"

// idiotic compile time fixed lens:
// static const char *lensfilename = "../optics/lenses/simple.fx";
// static const char *lensfilename = "../optics/lenses/CanonZoom.fx";
static const char *lensfilename = "../optics/lenses/tessar.fx";
static lens_element_t lenses[50];
static int lenses_cnt = 0;

#include <math.h>
#define CAMERA_BLADES 6

static inline void camera_display_info(camera_t *c)
{
  if(c->exposure_value > 6)
    display_print(rt.display, 0, 0, " %s   1/%.0f f/%.1f iso%d", lens_name, roundf(1.0/camera_exposure_time[c->exposure_value]),
        camera_f_stop[c->aperture_value+2], (int)c->iso);
  else
    display_print(rt.display, 0, 0, " %s   %.1f\" f/%.1f iso%d", lens_name, camera_exposure_time[c->exposure_value],
        camera_f_stop[c->aperture_value+2], (int)c->iso);
}

static inline void camera_print_info(FILE *f)
{
  fprintf(f, "camera   : ray traced lens `%s'\n", lensfilename);
  fprintf(f, "           %s\n", lens_extra);
  fprintf(f, "  pos    : %f %f %f\n", rt.cam->pos[0], rt.cam->pos[1], rt.cam->pos[2]);
  float dir[3] = {0, 0, 1.0};
  quaternion_transform(&(rt.cam->orient), dir);
  fprintf(f, "  dir    : %f %f %f\n", dir[0], dir[1], dir[2]);
  fprintf(f, "  quat   : %f, (%f %f %f)\n", rt.cam->orient.w, rt.cam->orient.x[0], rt.cam->orient.x[1], rt.cam->orient.x[2]);
  fprintf(f, "  focus  : %fmm (sensor offset)\n", rt.cam->focus);
  fprintf(f, "  crop   : %.1f\n", rt.cam->crop_factor);
  fprintf(f, "  film   : %dmm x %dmm\n", (int)(rt.cam->film_width*100.0), (int)(rt.cam->film_height*100.0 + 0.5f));
  fprintf(f, "  iso    : %d\n", (int)rt.cam->iso);
  if(rt.cam->exposure_value > 6)
    fprintf(f, "  Tv     : 1/%.0f\n", roundf(1.0/camera_exposure_time[rt.cam->exposure_value]));
  else
    fprintf(f, "  Tv     : %.1f''\n", camera_exposure_time[rt.cam->exposure_value]);
  fprintf(f, "  Av     : f/%.1f\n", camera_f_stop[rt.cam->aperture_value+2]);
  // maybe incorporate zoom into the polynomial at some point later.
  // fprintf(f, "  f      : %dmm\n", (int)(rt.cam->focal_length*100.0 + 0.5f));
}

static inline float camera_aperture_radius(const camera_t *c)
{
  return fminf(lens_aperture_housing_radius, lens_focal_length/(2.0f * camera_f_stop[c->aperture_value+2]));
}

static inline void camera_set_focus(camera_t *c, float dist)
{
  float mmdist = dist * 100.0f; // convert dm to mm
  // trace a couple of adjoint rays from there to the sensor and see where we need to put the sensor plane.
  float in[5], out[5];
  in[4] = out[4] = .5f; // wavelength 500nm
  float off = 0.0f;
  const int S = 5;
  for(int s=1;s<=S;s++)
  for(int k=0;k<2;k++)
  {
    // input starts on the scene facing pupil:
    in[k] = lens_outer_pupil_radius*.3f*s/(float)S; // don't go too far off-center
    in[1-k] = 0.0f;
    in[2+k] = in[k]/mmdist;
    in[3-k] = 0.0f;
    evaluate_reverse(lenses, lenses_cnt, 0.0f, in, out);
    off += out[k]/out[2+k];
    // fprintf(stderr, "in %f %f out %f %f\n", in[k], in[2+k], out[k], out[2+k]);
    // out[k+2] = -out[k+2];
    // lens_evaluate(out, in);
    // fprintf(stderr, "in %f %f out %f %f\n", in[k], in[2+k], out[k], out[2+k]);
    // fprintf(stderr, "estimate %d %d : %f\n", s, k, out[k]/out[2+k]);
  }
  off /= 2.0f*S; // average two guesses
  // the focus plane/sensor offset:
  fprintf(stderr, "camera focus shift  : %fmm (sensor offset)\n", off);
  // negative because of reverse direction
  if((off == off) && (off < 10.0f) && (off > -10.0f))
    c->focus = -off; // in mm
}

// copied from lens.h
static inline int lens_clip_aperture(const float x, const float y, const float radius, const int blades)
{ 
  // early out
  if(x*x + y*y > radius*radius) return 0;
  float xx = radius; 
  float yy = 0.0f;
  for(int b=1;b<blades+1;b++)
  {      
    float tmpx, tmpy;
    common_sincosf(2.0f*(float)M_PI/blades * b, &tmpy, &tmpx);
    tmpx *= radius;
    tmpy *= radius;
    const float normalx = xx + tmpx;
    const float normaly = yy + tmpy;
    float dot0 = (normalx)*(x-xx) + (normaly)*(y-yy);
    if(dot0 > 0.0f) return 0;
    xx = tmpx;
    yy = tmpy;
  }
  return 1;
}
static inline float lens_det_dxy_omega(const float *out)
{
  // the transform is: f(out[2], out[3]) = (out[2], out[3])/len (renormalise to l2 norm)
  // and len = 1./sqrtf(1+out[2]^2 + out[3]^2) so the determinant of that will be the square:
  return 1.0f/(1.0f + out[2]*out[2] + out[3]*out[3]);
}

static inline float camera_sample(const camera_t *c, path_t *p)
{
#if 0
  // hacky stats:
  static uint64_t samples = 0;
  static uint64_t samples_success = 0;
#pragma omp atomic
  samples++;
#endif
  // terrible hack, lack of initialization callback:
  if(lenses_cnt == 0)
  {
#pragma omp critical
    {
      if(lenses_cnt == 0)
      {
        lenses_cnt = system_configuration(lensfilename, lenses, 50);
        fprintf(stderr, "[lensrt] loading lens system with %d elements\n", lenses_cnt);
      }
    }
  }

  const float i = pointsampler(p, s_dim_image_x);
  const float j = pointsampler(p, s_dim_image_y);
  const float r1 = pointsampler(p, s_dim_aperture_x);
  const float r2 = pointsampler(p, s_dim_aperture_y);

  const float dm2mm = 100.0f;

  // [x,y,dx,dy,lambda] on the sensor, aperture, and outgoing pupil:
  float sen[5], ape[5], out[5];
  sen[4] = ape[4] = out[4] = p->lambda/1000.0f; // nanometers to micrometers
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);

  // sample point on the sensor in mm:
  sen[0] = dm2mm*(i-.5f)*c->film_width;
  sen[1] = dm2mm*(j-.5f)*c->film_height;
  // remember pixel coordinates:
  p->sensor.pixel_i = (sen[0]/(dm2mm*c->film_width)  + .5f)*c->width;
  p->sensor.pixel_j = (sen[1]/(dm2mm*c->film_height) + .5f)*c->height;

  // initial guess for the direction, pointing right at the aperture point:
  sen[2] = (2.0f*(r1-.5f)*lens_inner_pupil_radius - sen[0])/(lens_focal_length + c->focus);
  sen[3] = (2.0f*(r2-.5f)*lens_inner_pupil_radius - sen[1])/(lens_focal_length + c->focus);
  const float A = 4.0f*lens_inner_pupil_radius*lens_inner_pupil_radius;

  // move to beginning of optical system:
  sen[0] += sen[2] * c->focus;
  sen[1] += sen[3] * c->focus;

  const float zoom = 0.0f;
  int error = evaluate(lenses, lenses_cnt, zoom, sen, out);
  if(error) return 0.0f;
  // crop out by outgoing pupil and
  // crop at inward facing pupil:
  const float px = sen[0] + sen[2] * lens_focal_length, py = sen[1] + sen[3]*lens_focal_length;
  if((out[0]*out[0] + out[1]*out[1] > lens_outer_pupil_radius*lens_outer_pupil_radius) ||
     (px*px + py*py > lens_inner_pupil_radius*lens_inner_pupil_radius))
  {
    p->v[0].mode = s_absorb;
    return 0.0f;
  }

  // crop at aperture
  error = evaluate_aperture(lenses, lenses_cnt, zoom, sen, ape);
  if(error) return 0.0f;
  if(!lens_clip_aperture(ape[0], ape[1], aperture_radius, num_blades)) return 0.0f;


  // convert to world space:
  float camera_space_pos[3], camera_space_omega[3];
  camera_space_pos[0] = out[0]/dm2mm;
  camera_space_pos[1] = out[1]/dm2mm;
  camera_space_pos[2] = 0.0f;
  const float l2_norm = sqrtf(1.0f + out[2]*out[2] + out[3]*out[3]);
  camera_space_omega[0] = out[2]/l2_norm;
  camera_space_omega[1] = out[3]/l2_norm;
  camera_space_omega[2] = 1.0f/l2_norm;

  for(int k=0;k<3;k++)
  {
    p->v[0].hit.x[k] = c->pos[k] + camera_space_pos[0] * c->rg[k] + camera_space_pos[1] * c->up[k];
    p->e[1].omega[k] = camera_space_omega[0] * c->rg[k] + camera_space_omega[1] * c->up[k] + camera_space_omega[2] * c->dir[k];
  }

#if 0
#pragma omp atomic
  samples_success++;

  if(!(samples_success & 0xffff))
  {
    fprintf(stderr, "[lensrt] survival rate %.03f\n", (double)samples_success/(double)samples);
    samples_success = samples = 0;
  }
#endif

  // also clip to valid range.
  if(p->sensor.pixel_i < 0.0f || p->sensor.pixel_i > rt.width-1 ||
     p->sensor.pixel_j < 0.0f || p->sensor.pixel_j > rt.height-1)
  {
    p->v[0].mode = s_absorb;
    return 0.0f;
  }

  // compute geometric term and pdfs
  const float inv_p = c->film_width*c->film_height*A/(dm2mm*dm2mm);
  const float sensor = 100.0f*camera_exposure_time[c->exposure_value];

  // adjust pathspace parameters on vertices:
  p->v[0].hit.prim = -1;
  for(int k=0;k<3;k++)
    p->v[0].hit.n[k] = p->v[0].hit.gn[k] = c->dir[k];
  p->v[0].hit.shader = -1;
  p->v[0].pdf = 1.0f/inv_p;
  p->v[0].mode = s_sensor;
  p->v[0].rand_cnt = 6;

  p->v[1].rand_beg = p->v[0].rand_beg + p->v[0].rand_cnt;
  p->v[1].rand_cnt = 1; // only free path length
  p->v[1].pdf = 1.0f; // will be transformed by geometric term later on.

  // sensor/pdf * geometric term.
  return sensor*inv_p * lens_det_dxy_omega(sen) * lens_det_dxy_omega(sen)*dm2mm*dm2mm/(lens_focal_length*lens_focal_length);
}

static inline float camera_pdf(camera_t *c, path_t *p, int v)
{
  // TODO:
  return 0.0f;
}

static inline float camera_connect(const camera_t *c, path_t *p)
{
  // terrible hack, lack of initialization callback:
  if(lenses_cnt == 0)
  {
#pragma omp critical
    {
      if(lenses_cnt == 0)
      {
        lenses_cnt = system_configuration(lensfilename, lenses, 50);
        fprintf(stderr, "[lensrt] loading lens system with %d elements\n", lenses_cnt);
      }
    }
  }

  const int v = p->length; // create a new vertex at the end

  // construct lens hit:
  // const float ax = pointsampler(p, s_dim_light1);
  // const float ay = pointsampler(p, s_dim_light2);
  const float ox = pointsampler(p, s_dim_nee_x);
  const float oy = pointsampler(p, s_dim_nee_y);

  const float dm2mm = 100.0f;

  // [x,y,dx,dy,lambda] on the sensor, aperture, and outgoing pupil:
  float sen[5], ape[5], out[5];
  sen[4] = ape[4] = out[4] = p->lambda/1000.0f; // nanometers to micrometers

  // sample a point on the outgoing pupil:
  out[0] = cosf(2*M_PI*ox)*sqrtf(oy)*lens_outer_pupil_radius;
  out[1] = sinf(2*M_PI*ox)*sqrtf(oy)*lens_outer_pupil_radius;

  float view[3];
  for(int k=0;k<3;k++)
  {
    p->v[v].hit.x[k] = c->pos[k] + c->rg[k] * out[0] / dm2mm + c->up[k] * out[1] / dm2mm;
    p->v[v].hit.a[k] = c->rg[k];
    p->v[v].hit.b[k] = c->up[k];
    p->v[v].hit.gn[k] = p->v[v].hit.n[k] = c->dir[k];
    view[k] = p->v[v-1].hit.x[k] - p->v[v].hit.x[k];
  }
  p->e[v].dist = sqrtf(dotproduct(view, view));
  for(int k=0;k<3;k++) p->e[v].omega[k] = - view[k] / p->e[v].dist;
  // project to two plane light field parametrisation:
  const float norm = -dotproduct(c->dir, view);
  out[2] = dotproduct(c->rg, view)/norm;
  out[3] = dotproduct(c->up, view)/norm;

  // compute outgoing direction:
  const float zoom = 0.0f;
  int error = evaluate_reverse(lenses, lenses_cnt, zoom, out, sen);
  // fprintf(stderr, "eval: %f %f %f %f -- %f %f %f %f\n", out[0], out[1], out[2], out[3], sen[0], sen[1], sen[2], sen[3]);
  if(error) return 0.0f;

  // swap the sign for aperture computation (reverse path now)
  // we depend on that later when computing pixel_{i,j}, too.
  sen[2] = -sen[2];
  sen[3] = -sen[3];
#if 1
  error = evaluate_aperture(lenses, lenses_cnt, zoom, sen, ape);
  if(error) return 0.0f;
#endif
  const int num_blades = CAMERA_BLADES;
  const float aperture_radius = camera_aperture_radius(c);
  // const float A = lens_aperture_area(aperture_radius, num_blades);

  // minus sensor direction because we swapeed the sign for aperture computation
  p->sensor.pixel_i = .5f*c->width  + (sen[0]-sen[2]*c->focus)/(dm2mm*c->film_width) * c->width;
  p->sensor.pixel_j = .5f*c->height + (sen[1]-sen[3]*c->focus)/(dm2mm*c->film_height) * c->height;
  const float sensor = 100.0f*camera_exposure_time[c->exposure_value];
  p->v[v].pdf = dm2mm*dm2mm/(M_PI*lens_outer_pupil_radius*lens_outer_pupil_radius);
  p->v[v].mode = s_sensor;
  p->v[v].flags = s_none;
  p->v[v].hit.prim = -1;
  p->v[v].hit.shader = -1;
  p->v[v].shading.roughness = 1.0f; // instruct path space that we don't need an additional cosine

  if(!lens_clip_aperture(ape[0], ape[1], aperture_radius, num_blades)) return 0.0f;
#if 1
  // already done in ray tracing:
  // crop at inward facing pupil:
  // const float px = sen[0] + sen[2] * lens_focal_length, py = sen[1] + sen[3]*lens_focal_length;
  // if(px*px + py*py > lens_inner_pupil_radius*lens_inner_pupil_radius) return 0.0f;
#endif
  if(!(p->pixel_i >= 0 && p->pixel_j >= 0 && p->pixel_i < c->width && p->pixel_j < c->height)) return 0.0f;
  return sensor/(p->v[v].pdf);
}

static inline float camera_pdf_connect(const camera_t *c, path_t *p, int v)
{
  // sampled the outward facing pupil:
  const float dm2mm = 100.0f;
  return dm2mm*dm2mm/(M_PI*lens_outer_pupil_radius*lens_outer_pupil_radius);
}

static inline float camera_eval(const camera_t *c, path_t *p)
{
  return 1.0f;
}


#endif
