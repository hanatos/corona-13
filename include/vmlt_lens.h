#ifndef VMLT_LENS_H
#define VMLT_LENS_H

#include "points.h"
#include "pathspace/manifold.h"

void *lens_init()
{
  return 0;
}

void lens_cleanup(void *data) { }

float lens_suitability(const path_t *p, void *data)
{
  // currently only implemented for paths started at the eye
  if(!(p->v[0].mode & s_sensor)) return 0.0f;
  // only for purely specular paths, to explore highlights
  for(int k=1;k<p->length-1;k++) if(!(p->v[k].mode & s_specular)) return 0.0f;
  return 1.0f;
}

static int lens_reflect(path_t *p, int v)
{
  float H[3]; // reconstruct world space half vector
  // half vector lives in ortho normal basis hit.{a,b,n}
  for(int k=0;k<3;k++)
    H[k] = p->v[v].hit.a[k] * p->v[v].diffgeo.h[0] +
           p->v[v].hit.b[k] * p->v[v].diffgeo.h[1] +
           p->v[v].hit.n[k] * p->v[v].diffgeo.h[2];

  // now update outgoing direction accordingly
  if(p->v[v].mode & s_reflect)
  {
    const float dot = dotproduct(p->e[v].omega, H);
    for(int k=0;k<3;k++)
      p->e[v+1].omega[k] = p->e[v].omega[k] - 2.0f*H[k] * dot;
  }
  else if(p->v[v].mode & s_transmit)
  {
    const float eta_ratio = path_eta_ratio(p, v); // = n1/n2;
    if(fabsf(eta_ratio - 1.0f) > 1e-5f)
    {
      float cosr = dotproduct(H, p->e[v].omega);
      if(cosr > 0.0f) // need to flip H to point into R hemisphere
        for(int k=0;k<3;k++) H[k] = -H[k];
      else cosr = -cosr;
      if(eta_ratio < 0.0f)
        return 1; // volume nesting broken.
      const float cost2 = 1.0f - eta_ratio*eta_ratio * (1.0f - cosr*cosr);
      if(cost2 <= 0.0f)
        return 1; // total internal reflection
      const float cost = sqrtf(cost2);
      const float f = eta_ratio*cosr - cost;
      for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k]*eta_ratio + f * H[k];
      normalise(p->e[v+1].omega);
    }
    else
    {
      for(int k=0;k<3;k++) p->e[v+1].omega[k] = p->e[v].omega[k];
    }
  }
  else if(p->v[v].mode & s_volume)
  {
    for(int k=0;k<3;k++)
      p->e[v+1].omega[k] = H[k];
  }
  return 0;
}

float lens_mutate(path_t *curr, path_t *tent, void *data)
{
  const int n = curr->length-1;

  curr->v[0].diffgeo.type = s_pinned_position;
  for(int i=1;i<=n;i++)
    curr->v[i].diffgeo.type = s_free;
  manifold_compute_tangents(curr, 0, n); // only reflect_H above depends on that..
  memcpy(tent, curr, sizeof(path_t)); // wasteful, but whatever.

  // replace that by _thread at some point
  const int tid = common_get_threadid();

  // TODO: also mutate wavelength a bit?
  // tent->lambda = spectrum_mutate(curr->lambda, points_rand(rt.points, tid), 0);

  float pdf;
  // trace tent path keeping half vectors to re-initialise
  // all hit point infos, medium transactions and to get new error vector

  // mutate first direction
#if 1 // stupid hg mutation
  float a[3], b[3];
  float *w = tent->e[1].omega;
  const float g = 0.995f;
  const float sqr = (1.0f-g*g)/(1.0f-g+2.0f*g*points_rand(rt.points, tid));
  const float cos_theta = 1.0f/(2.0f*g) * (1.0f + g*g - sqr*sqr);
  const float phi = 2.0f*M_PI*points_rand(rt.points, tid);
  const float l = sqrtf(1.0f-cos_theta*cos_theta);
  const float aa = cosf(phi)*l, bb = sinf(phi)*l;
  get_onb(w, a, b);
  for(int k=0;k<3;k++)
    w[k] = curr->e[1].omega[k] * cos_theta + a[k] * aa + b[k] * bb;
#endif
  for(int k=0;k<3;k++)
    tent->v[1].hit.x[k] = tent->v[0].hit.x[k] + tent->e[1].dist * tent->e[1].omega[k];
  if(path_project(tent, 1, &pdf) ||
      (tent->v[1].flags  != curr->v[1].flags) ||
      (primid_invalid(tent->v[1].hit.prim) != primid_invalid(curr->v[1].hit.prim))) 
    return 0.0f;

  // the others keep the half vectors constant
  for(int v=2;v<=n;v++)
  {
    if(lens_reflect(tent, v-1))
      return 0.0f;

    if(path_propagate(tent, v, &pdf) ||
        (tent->v[v-1].mode != curr->v[v-1].mode) ||
        (tent->v[v].flags  != curr->v[v].flags) ||
        (primid_invalid(tent->v[v].hit.prim) != primid_invalid(curr->v[v].hit.prim))) 
      return 0.0f;
  }

  const double f_tent = path_measurement_contribution_dwp(tent);
  const double f_curr = path_measurement_contribution_dwp(curr);
  return f_tent/f_curr;
}

void lens_print_info(FILE *f, void *data)
{
  fprintf(f, "         : lens perturbation\n");
}

#endif
