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

#include "camera.h"
#include "pathspace.h"
#include "pointsampler.h"
#include "spectrum.h"
#include "view.h"
#include "display.h"

#include <math.h>

// convert X+Y+Z=1 input light source to something visible:
static const float camera_sensor_response = 106.86535f;

int camera_global_init() { return 0; }

void camera_display_info(const camera_t *c)
{
  if(c->exposure_value > 6)
    display_print(rt.display, 0, 0, " 1/%.0f f/%.1f %.0fmm iso%d", roundf(1.0/view_tv2time(c->exposure_value)),
        view_av2fstop(c->aperture_value), 100.0*c->focal_length, (int)c->iso);
  else
    display_print(rt.display, 0, 0, " %.1f\" f/%.1f %.0fmm iso%d", view_tv2time(c->exposure_value),
        view_av2fstop(c->aperture_value), 100.0*c->focal_length, (int)c->iso);
}

void camera_print_info(const camera_t *c, FILE *f)
{
  fprintf(f, "camera   : thin lens model\n");
  fprintf(f, "  focus  : %f\n", c->focus);
  fprintf(f, "  film   : %dmm x %dmm\n", (int)(c->film_width*100.0f+0.5f), (int)(c->film_height*100.0f + 0.5f));
  if(c->exposure_value > 6)
    fprintf(f, "         : 1/%.0f f/%.1f %.0fmm iso %d\n", roundf(1.0/view_tv2time(c->exposure_value)),
        view_av2fstop(c->aperture_value), 100.0*c->focal_length, (int)c->iso);
  else
    fprintf(f, "         : %.1f\" f/%.1f %.0fmm iso %d\n", view_tv2time(c->exposure_value),
        view_av2fstop(c->aperture_value), 100.0*c->focal_length, (int)c->iso);
}

void camera_set_focus(camera_t *c, float dist)
{
  c->focus = dist;
}

static float _camera_aperture_area(const camera_t *c)
{
  const float n = view_av2fstop(c->aperture_value);
  const float f = c->focal_length;
  const float A = M_PI*f*f/(4.0f*n*n);
  return A;
}

static mf_t _camera_sample_internal(
    const camera_t *c,
    path_t *p,
    const float i,
    const float j,
    const float u,
    const float v)
{
  // init tangent frame
  view_cam_init_frame(p, &p->v[0].hit);

  const float f = c->focus / c->focal_length;
  const float f_dir =   c->focus;
  const float f_rg  = - c->film_width  * f / view_width();
  const float f_up  = - c->film_height * f / view_height();

  float aoff[3];
  for(int k=0;k<3;k++) aoff[k] = u * p->v[0].hit.a[k] + v * p->v[0].hit.b[k];
  for(int k=0;k<3;k++) p->e[1].omega[k] = f_dir*p->v[0].hit.n[k] + ((i-.5f*view_width())*f_rg*p->v[0].hit.a[k] + (j-.5f*view_height())*f_up*p->v[0].hit.b[k]) - aoff[k];
  normalise(p->e[1].omega);

  const float A = _camera_aperture_area(c);
  const float pdf_a = 1./A; // vertex area pdf of vertex[0] on the aperture
  const float sensor = camera_sensor_response*100.0f*view_tv2time(c->exposure_value);
  const float dot = dotproduct(p->e[1].omega, p->v[0].hit.n); // cos between aperture and focal plane
  const float dot4 = dot*dot*dot*dot;
  p->v[0].hit.prim = INVALID_PRIMID;
  p->v[0].hit.shader = -1;
  p->v[0].pdf = mf_set1(1.0f); // first two vertices created in one go, pdf is only in v[1]
  p->v[0].material_modes = p->v[0].mode = s_sensor;
  p->sensor.pixel_i = CLAMP(i, 0.0, view_width()-1e-4f);
  p->sensor.pixel_j = CLAMP(j, 0.0, view_height()-1e-4f);
  p->v[0].rand_cnt = s_dim_num_pt_beg;

  p->v[1].rand_beg = p->v[0].rand_beg + p->v[0].rand_cnt;
  p->v[1].rand_cnt = 1; // only free path length
  const float G = dot4/(c->focal_length * c->focal_length);  // geometric term up to point on the plane in focus (outside lens), focus distance will cancel out later
  // pdf of outgoing ray in projected solid angle measure = p_v / G
  const float pdf_v = 1.0f/(c->film_width*c->film_height); // surface area pdf on film
  p->v[1].pdf = mf_set1(pdf_v * pdf_a / G); // will be transformed to vertex area in scene by geometric term later on.
  // assert(fabsf(camera_pdf(c, p, 0) - p->v[1].pdf) < 0.001f*fabsf(p->v[1].pdf));

  for(int k=0;k<3;k++) p->v[0].hit.x[k] += aoff[k];
  // throughput is measurement contribution (only geometric term here) times sensor response divided by 1/sensor and 1/aperture
  return mf_set1(sensor * G / (pdf_a * pdf_v));
}

mf_t camera_sample(const camera_t *c, path_t *p)
{
  const float i = p->sensor.pixel_set ? p->sensor.pixel_i : (pointsampler(p, s_dim_image_x)*view_width());
  const float j = p->sensor.pixel_set ? p->sensor.pixel_j : (pointsampler(p, s_dim_image_y)*view_height());

  const float r1 = p->sensor.aperture_set ? p->sensor.aperture_x : pointsampler(p, s_dim_aperture_x);
  const float r2 = p->sensor.aperture_set ? p->sensor.aperture_y : pointsampler(p, s_dim_aperture_y);
  p->sensor.aperture_x = r1;
  p->sensor.aperture_y = r2;
  const float lens_radius = (.5f/view_av2fstop(c->aperture_value))*c->focal_length;
  float u = cosf(2*M_PI*r1)*sqrtf(r2)*lens_radius;
  float v = sinf(2*M_PI*r1)*sqrtf(r2)*lens_radius;
  return _camera_sample_internal(c, p, i, j, u, v);
}

mf_t camera_mutate_aperture(const camera_t *c, path_t *p, const float r1, const float r2, const float step)
{
  // reconstruct current vertex coordinates on lens:
  const float q1 = p->sensor.aperture_x;
  const float q2 = p->sensor.aperture_y;
  const float lens_radius = (.5f/view_av2fstop(c->aperture_value))*c->focal_length;
  const float u = cosf(2*M_PI*q1)*sqrtf(q2)*lens_radius;
  const float v = sinf(2*M_PI*q1)*sqrtf(q2)*lens_radius;
  // mutate:
  const float uu = (.5f - r1) * lens_radius * step;
  const float vv = (.5f - r2) * lens_radius * step;
  return _camera_sample_internal(c, p, p->sensor.pixel_i, p->sensor.pixel_j, u + uu, v + vv);
}

// vertex area pdf(tent->v[0] | curr->v[0])
mf_t camera_pdf_mutate_aperture(const camera_t *c, const path_t *curr, const path_t *tent, const float step)
{
  float x0c[3], x0t[3];
  for(int k=0;k<3;k++) x0c[k] = c->pos[k] * (1.0f-curr->time) + c->pos_t1[k] * curr->time;
  for(int k=0;k<3;k++) x0t[k] = c->pos[k] * (1.0f-tent->time) + c->pos_t1[k] * tent->time;
  for(int k=0;k<3;k++) x0c[k] -= curr->v[0].hit.x[k];
  for(int k=0;k<3;k++) x0t[k] -= tent->v[0].hit.x[k];
  const float uc = dotproduct(x0c, curr->v[0].hit.a);
  const float vc = dotproduct(x0c, curr->v[0].hit.b);
  const float ut = dotproduct(x0t, tent->v[0].hit.a);
  const float vt = dotproduct(x0t, tent->v[0].hit.b);
  const float lens_radius = (.5f/view_av2fstop(c->aperture_value))*c->focal_length;
  if(fabsf(uc - ut) > .5 * step * lens_radius) return mf_set1(0.0f);
  if(fabsf(vc - vt) > .5 * step * lens_radius) return mf_set1(0.0f);
  return mf_set1(1.0f/(step*step * lens_radius*lens_radius));
}

// belongs to sample(), returns pdf in projected outgoing solid angle of p->e[{0,v}].omega (depending on trace direction)
mf_t camera_pdf(const camera_t *c, const path_t *p, int v)
{
  // either produce v[0] and v[1] or would have produced v[len-1] and v[len-2]
  assert(v == 0 || v == p->length-1);
  // includes both p(x) and p(y) on sensor and aperture
  const float A = _camera_aperture_area(c);
  const float pdf_v = 1./(A*c->film_width*c->film_height);
  // now transform to pdf_omega by dividing out the geometric term:
  // use direction of edge leading up to vertex v from outside the camera.
  const int e = v ? v : 1;
  const float dot = fabsf(dotproduct(p->e[e].omega, p->v[v].hit.n));
  const float G = (dot*dot*dot*dot)/(c->focal_length * c->focal_length);
  // pdf of outgoing ray in projected solid angle measure = p_v / G = 1/A * cos^4/f^2
  // transform to projected solid angle measure, will be transformed to on-surface probability in pathspace.
  return mf_set1(pdf_v / G);
}

mf_t camera_connect(const camera_t *c, path_t *p)
{
  //             path tracing case               mutation case, end vertex there already              light tracing case, create a new vertex at the end     
  const int v0 = (p->v[0].mode & s_sensor) ? 0 : ((p->v[p->length-1].mode & s_sensor) ? p->length-1 : p->length);
  const int v1 = (p->v[0].mode & s_sensor) ? 1 : ((p->v[p->length-1].mode & s_sensor) ? p->length-2 : p->length-1);
  const int e  = (p->v[0].mode & s_sensor) ? 1 : ((p->v[p->length-1].mode & s_sensor) ? p->length-1 : p->length);
  // construct lens hit:
  const float x = p->sensor.aperture_set ? p->sensor.aperture_x : pointsampler(p, s_dim_nee_x);
  const float y = p->sensor.aperture_set ? p->sensor.aperture_y : pointsampler(p, s_dim_nee_y);
  float offs[3], view[3];
  const float lens_radius = (.5f/view_av2fstop(c->aperture_value))*c->focal_length;
  float r = cosf(2*M_PI*x)*sqrtf(y)*lens_radius;
  float s = sinf(2*M_PI*x)*sqrtf(y)*lens_radius;
  view_cam_init_frame(p, &p->v[v0].hit);
  for(int k=0;k<3;k++)
  {
    offs[k] = r * p->v[v0].hit.a[k] + s * p->v[v0].hit.b[k];
    p->v[v0].hit.x[k] += offs[k];
    view[k] = p->v[v1].hit.x[k] - p->v[v0].hit.x[k];
  }
  if(dotproduct(p->v[v0].hit.gn, view) <= 0) return mf_set1(0.0f);
  p->e[e].dist = sqrtf(dotproduct(view, view));
  if(p->v[0].mode & s_sensor)
    for(int k=0;k<3;k++) p->e[e].omega[k] = view[k] / p->e[e].dist;
  else
    for(int k=0;k<3;k++) p->e[e].omega[k] = - view[k] / p->e[e].dist;

  // project onto basis vectors: first transform v to focal plane along ray to sampled lens point.
  // then let v be vector to center of lens, subtract basis vector in dir and project to rg, up.
  float dot = dotproduct(p->v[v0].hit.n, view);
  const float f = c->focus / c->focal_length;
  const float f_dir =   c->focus;
  const float f2 = f_dir/dot;
  const float f_rg  = - c->film_width  * f / view_width();
  const float f_up  = - c->film_height * f / view_height();
  for(int k=0;k<3;k++) view[k] = view[k]*f2 + offs[k] - p->v[v0].hit.n[k]*f_dir;
  p->sensor.pixel_i = .5f*view_width()  + dotproduct(view, p->v[v0].hit.a)/f_rg;
  p->sensor.pixel_j = .5f*view_height() + dotproduct(view, p->v[v0].hit.b)/f_up;
  p->sensor.aperture_x = x;
  p->sensor.aperture_y = y;
  const float A = _camera_aperture_area(c);
  const float sensor = camera_sensor_response*100.0f*view_tv2time(c->exposure_value);
  p->v[v0].pdf = mf_set1(1./A);
  p->v[v0].material_modes =
  p->v[v0].mode = s_sensor;
  p->v[v0].flags = s_none;
  p->v[v0].hit.prim = INVALID_PRIMID;
  p->v[v0].hit.shader = -1;
  p->v[v0].shading.roughness = 1.0f; // instruct path space that we don't need an additional cosine

  // cos from position on pixel plane to center of aperture
  // (using the cos from position on pixel plane (omega.x -= r, omega.y -= s)
  // to actual position on aperture would result in a slightly too dark picture)
  float omega0[3] = { (p->sensor.pixel_i - 0.5f*view_width()) / view_width() * c->film_width,
                      (p->sensor.pixel_j - 0.5f*view_height()) / view_height() * c->film_height,
                       c->focal_length};
  normalise(omega0);

  if(!(p->sensor.pixel_i >= 0 && p->sensor.pixel_j >= 0 && p->sensor.pixel_i < view_width() && p->sensor.pixel_j < view_height())) return mf_set1(0.0f);
  // throughput is sensor response / pdf_a (p(x) cancels out with the G term)
  return mf_set1(sensor * A);
}

mf_t camera_pdf_connect(const camera_t *c, const path_t *p, int v)
{
  // sampled only the aperture surface, return 1/A * p(x) on sensor
  const float A = _camera_aperture_area(c);
  return mf_set1(1./A);
}

mf_t camera_eval(const camera_t *c, const path_t *p)
{
  if(!(p->sensor.pixel_i >= 0 && p->sensor.pixel_j >= 0 && p->sensor.pixel_i < view_width() && p->sensor.pixel_j < view_height())) return mf_set1(0.0f);
  const float sensor = camera_sensor_response*100.0f*view_tv2time(c->exposure_value);
  return mf_set1(sensor);
}

// worst technique of all, never used so far.
#if 0
static inline mf_t camera_intersect(const camera_t *c, const ray_t *ray, const rayhit_t *hit, float *i, float *j, float *t)
{
  const float lambda2 = hit->lambda*hit->lambda;
  const float lambda0_2 = 525.0f*525.0f;
  const float lca = ((lambda0_2 + c->chromatic_aberration)/(lambda0_2)) * lambda2/(lambda2 + c->chromatic_aberration);
  const float tca = (lca + c->transverse_chromatic_aberration)/(c->transverse_chromatic_aberration + 1.0f);
  const float dot = dotproduct(ray->dir, c->dir);
  if(dot >= 0.0f) { *t = INFINITY; return 0.0f; }
  float v[3], offs[3];
  // intersect outgoing pupil with ray, get point on lens
  for(int k=0;k<3;k++) v[k] = c->pos[k] - ray->pos[k];
  *t = dotproduct(v, c->dir)/dotproduct(ray->dir, c->dir);
  if(*t < 0) { *t = INFINITY; return 0.0f; }
  for(int k=0;k<3;k++) offs[k] = ray->pos[k] + *t*ray->dir[k] - c->pos[k];
  if(dotproduct(offs, offs) > c->lens_radius*c->lens_radius) { *t = INFINITY; return 0.0f; }
  // project onto basis vectors.
  // const float f = c->focus/dotproduct(c->dir, ray->dir);
  // for(int k=0;k<3;k++) v[k] = ray->dir[k]*f + offs[k] - c->f_ul[k];
  // const float inv_fr_len = (c->width*c->focal_length)/(c->film_width * c->focus);
  // const float inv_fd_len = (c->height*c->focal_length)/(c->film_height * c->focus);
  // *i = dotproduct(v, c->f_r)*inv_fr_len*inv_fr_len;
  // *j = dotproduct(v, c->f_d)*inv_fd_len*inv_fd_len;
  const float f = lca*c->f_dir/dotproduct(c->dir, v);
  for(int k=0;k<3;k++) v[k] = v[k]*f + offs[k] - c->dir[k]*lca*c->f_dir;
  *i = .5f*c->width  + dotproduct(v, c->rg)/(tca*c->f_rg);
  *j = .5f*c->height + dotproduct(v, c->up)/(tca*c->f_up);
  if(*i < 0 || *j < 0 || *i >= c->width || *j >= c->height) { *t = INFINITY; return 0.0f; }
  const float n = camera_f_stop[c->aperture_value];
  const float sensor = camera_sensor_response*100.0f*camera_exposure_time[c->exposure_value]*M_PI/(4.0f*n*n);
  return sensor;
}
#endif
