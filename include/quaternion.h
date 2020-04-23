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

// copyright (c) 2007 j.hanika
// GPL

#pragma once
#include <math.h>

typedef struct
{
  float w, x[3];
}
quaternion_t;

static inline void quaternion_init(quaternion_t *qt, const float angle, const float * const axis)
{
  float tmp = sinf(angle/2.0f);
  qt->w = cosf(angle/2.0f);
  for(int k=0;k<3;k++) qt->x[k] = tmp * axis[k];
}

static inline void quaternion_invert(quaternion_t *q)
{
  for(int k=0;k<3;k++) q->x[k] = -q->x[k];
}

static inline void quaternion_mult(quaternion_t *in, const quaternion_t *p)
{
  quaternion_t r = *in;
  in->x[0] = r.w*p->x[0] + r.x[0]*p->w    + r.x[1]*p->x[2] - r.x[2]*p->x[1];
  in->x[1] = r.w*p->x[1] - r.x[0]*p->x[2] + r.x[1]*p->w    + r.x[2]*p->x[0];
  in->x[2] = r.w*p->x[2] + r.x[0]*p->x[1] - r.x[1]*p->x[0] + r.x[2]*p->w;
  in->w    = r.w*p->w    - r.x[0]*p->x[0] - r.x[1]*p->x[1] - r.x[2]*p->x[2];
}

static inline void quaternion_mult_fleft(const quaternion_t *r, quaternion_t *in)
{
  quaternion_t p = *in;
  in->x[0] = r->w*p.x[0] + r->x[0]*p.w    + r->x[1]*p.x[2] - r->x[2]*p.x[1];
  in->x[1] = r->w*p.x[1] - r->x[0]*p.x[2] + r->x[1]*p.w    + r->x[2]*p.x[0];
  in->x[2] = r->w*p.x[2] + r->x[0]*p.x[1] - r->x[1]*p.x[0] + r->x[2]*p.w;
  in->w    = r->w*p.w    - r->x[0]*p.x[0] - r->x[1]*p.x[1] - r->x[2]*p.x[2];
}

static inline void quaternion_transform(const quaternion_t *const q, float* p)
{
  // q * v * q'
  quaternion_t vectorQuat, inverseQuat, resultQuat;

  for(int k=0;k<3;k++) vectorQuat.x[k] = p[k];
  vectorQuat.w = 0.0f;

  inverseQuat = *q;
  quaternion_invert(&inverseQuat);

  resultQuat = *q;
  quaternion_mult(&resultQuat, &vectorQuat);
  quaternion_mult(&resultQuat, &inverseQuat);
  for(int k=0;k<3;k++)
    p[k] = resultQuat.x[k];
}

static inline void quaternion_lerp(const quaternion_t *const q, const quaternion_t *const p, const float t, quaternion_t *res)
{
  res->w = (1.0f-t)*q->w + t*p->w;
  for(int k=0;k<3;k++)
    res->x[k] = (1.0f-t)*q->x[k] + t*p->x[k];
  // TODO: i think the full quat should be re-normalised here.
}

static inline void quaternion_slerp(const quaternion_t *const q, const quaternion_t *const p, const float t, quaternion_t *res)
{
  // explicitly expanded dotproduct(q->x, p->x) so we don't depend on any other
  // headers.  this is a good idea because external tools can then write
  // cameras which composite quaternions without loading too much stuff.
  float cos_theta_2 = q->w * p->w + (q->x[0]*p->x[0] + q->x[1]*p->x[1] + q->x[2]*p->x[2]);
  if(fabsf(cos_theta_2) >= 1.0f)
  {
    *res = *q;
    return;
  }
  float theta_2 = acosf(cos_theta_2);
  float sin_theta_2 = sqrtf(1.0f - cos_theta_2*cos_theta_2);
  if(fabsf(sin_theta_2) < 1e-10f)
  {
    res->w = (q->w + p->w)*.5f;
    for(int k=0;k<3;k++) res->x[k] = (q->x[k] + p->x[k])*.5f;
    return;
  }
  float a = sinf((1.0f - t) * theta_2) / sin_theta_2;
  float b = sinf(t * theta_2) / sin_theta_2;
  res->w = q->w * a + p->w * b;
  for(int k=0;k<3;k++)
    res->x[k] = q->x[k] * a + p->x[k] * b;
}
