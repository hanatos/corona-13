#ifndef CORONA_GEO_SHELL_PROC_H
#define CORONA_GEO_SHELL_PROC_H

// simple analytic sphere for debugging

#if 0 // sphere:
static inline float _geo_shell_tex_min_free_path(const prims_t *p, const primid_t pi, const float s, const float t, const float w)
{
  const float sx = 1.f, sy = 1.0f;
  const float x = (sx*s - floorf(sx*s));
  const float y = (sy*t - floorf(sy*t));
  const float z = w;
  const float cx = .5f;
  const float cy = .5f;
  const float cz = .5f;
  const float rad = 0.3f;
  return sqrtf((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz)) - rad;
}

static inline void _geo_shell_tex_get_normal(const prims_t *p, const primid_t pi, const float s, const float t, const float w, float *n)
{
  const float sx = 1.f, sy = 1.0f;
  const float x = (sx*s - floorf(sx*s));
  const float y = (sy*t - floorf(sy*t));
  const float z = w;
  const float cx = .5f;
  const float cy = .5f;
  const float cz = .5f;
  // const float rad = 0.3f;
  n[0] = .5f/sqrtf((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz)) * 2.0f*(x-cx);
  n[1] = .5f/sqrtf((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz)) * 2.0f*(y-cy);
  n[2] = .5f/sqrtf((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz)) * 2.0f*(z-cz);
}
#endif

#if 1 // chain armour
static inline float _geo_shell_tex_min_free_path(const prims_t *p, const primid_t pi, const float s, const float t, const float w)
{
  // two torii defined on the torus :)
  const float c0[3] = {0.0f, 0.0f, 0.5f};
  const float c1[3] = {0.5f, 0.5f, 0.5f};
  const float r = 0.05f, R = 0.42f;
  const float tilt = 3.0f;
  const float nn = sqrtf(tilt*tilt*r*r + R*R);
  // TODO: perturb normals slightly based on absolute texture coords :)
  const float n0[3] = {0.0f, tilt*r/nn, R/nn};
  const float n1[3] = {0.0f, -tilt*r/nn, R/nn};

  const float sx = 1.0f, sy = 1.0f;
  float tx, ty, tz = w;
  float dist = FLT_MAX;
  for(int k=0;k<4;k++)
  {
    tx = (sx*s - floorf(sx*s));
    ty = (sy*t - floorf(sy*t));
    if(k&1) tx -= 1.0f;
    if(k&2) ty -= 1.0f;

    const float xt0[3] = {tx-c0[0], ty-c0[1], tz-c0[2]};
    float a0[3], b0[3];
    get_onb(n0, a0, b0);
    const float x0[3] = {dotproduct(xt0, a0), dotproduct(xt0, b0), dotproduct(xt0, n0)};
    const float norm0 = sqrtf(x0[0]*x0[0] + x0[1]*x0[1]);
    const float dd0[3] = {x0[0]/norm0 * R - x0[0], x0[1]/norm0 * R - x0[1], 0.0f - x0[2]};
    dist = fminf(dist, sqrtf(dotproduct(dd0, dd0)) - r);

    const float xt1[3] = {tx-c1[0], ty-c1[1], tz-c1[2]};
    float a1[3], b1[3];
    get_onb(n1, a1, b1);
    const float x1[3] = {dotproduct(xt1, a1), dotproduct(xt1, b1), dotproduct(xt1, n1)};
    const float norm1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]);
    const float dd1[3] = {x1[0]/norm1 * R - x1[0], x1[1]/norm1 * R - x1[1], 0.0f - x1[2]};
    dist = fminf(dist, sqrtf(dotproduct(dd1, dd1)) - r);
  }
  return dist;
}

static inline void _geo_shell_tex_get_normal(const prims_t *p, const primid_t pi, const float s, const float t, const float w, float *n)
{
  memset(n, 0, sizeof(float)*3);
  const float c0[3] = {0.0f, 0.0f, 0.5f};
  const float c1[3] = {0.5f, 0.5f, 0.5f};
  const float r = 0.05f, R = 0.42f;
  const float tilt = 3.0f;
  const float nn = sqrtf(tilt*tilt*r*r + R*R);
  // TODO: perturb normals slightly based on absolute texture coords :)
  const float n0[3] = {0.0f, tilt*r/nn, R/nn};
  const float n1[3] = {0.0f, -tilt*r/nn, R/nn};

  const float sx = 1.0f, sy = 1.0f;
  float tx, ty, tz = w;
  float dist = FLT_MAX;
  float mdist;
  for(int k=0;k<4;k++)
  {
    tx = (sx*s - floorf(sx*s));
    ty = (sy*t - floorf(sy*t));
    if(k&1) tx -= 1.0f;
    if(k&2) ty -= 1.0f;

    const float xt0[3] = {tx-c0[0], ty-c0[1], tz-c0[2]};
    float a0[3], b0[3];
    get_onb(n0, a0, b0);
    const float x0[3] = {dotproduct(xt0, a0), dotproduct(xt0, b0), dotproduct(xt0, n0)};
    const float norm0 = sqrtf(x0[0]*x0[0] + x0[1]*x0[1]);
    const float dd0[3] = {x0[0]/norm0 * R - x0[0], x0[1]/norm0 * R - x0[1], 0.0f - x0[2]};
    mdist = sqrtf(dotproduct(dd0, dd0)) - r;
    if(mdist < dist)
    {
      for(int i=0;i<3;i++) n[i] = dd0[i]/(mdist + r);
      dist = mdist;
    }

    const float xt1[3] = {tx-c1[0], ty-c1[1], tz-c1[2]};
    float a1[3], b1[3];
    get_onb(n1, a1, b1);
    const float x1[3] = {dotproduct(xt1, a1), dotproduct(xt1, b1), dotproduct(xt1, n1)};
    const float norm1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]);
    const float dd1[3] = {x1[0]/norm1 * R - x1[0], x1[1]/norm1 * R - x1[1], 0.0f - x1[2]};
    mdist = sqrtf(dotproduct(dd1, dd1)) - r;
    if(mdist < dist)
    {
      for(int i=0;i<3;i++) n[i] = dd1[i]/(mdist + r);
      dist = mdist;
    }
  }
}

#endif

#endif
