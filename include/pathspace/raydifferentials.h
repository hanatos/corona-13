#ifndef _CORONA_PATHSPACE_RAYDIFFERENTIALS_H
#define _CORONA_PATHSPACE_RAYDIFFERENTIALS_H

#include "matrix2d.h"
#include "view.h"

static inline float _rd_inv(const path_t *p, const int v, const float *m, float *inv)
{ // wrapper for fast path for 2x2 blocks
  if(p->v[v].mode & s_volume) return mat3_invert(m, inv);
  else return mat3_invert_sub2(m, inv);
}

// returns the offset vectors in world space of v[1] relative to vertex v[1]
// of the given path which result if pixel (i+stepi,j) and (i,j+stepj) are
// sampled.
// returns != 0 in fail cases (path vertex is outside camera's viewing frustum)
static inline int raydifferentials_v1(
    const path_t *path,
    const float stepi,
    const float stepj,
    float *rd_i,
    float *rd_j)
{
  assert(path->v[0].mode & s_sensor);
  if(path->sensor.pixel_i < 0.0 || path->sensor.pixel_i >= view_width() ||
     path->sensor.pixel_j < 0.0 || path->sensor.pixel_j >= view_height())
    return 1;
  path_t p;
  float o[3], dotpo, dotpd, t;
  p.time = path->time;
  p.lambda = path->lambda;
  p.sensor = path->sensor;
  p.sensor.pixel_set = 1;
  p.sensor.aperture_set = 1;
  p.sensor.camid_set = 1;
  // sample new pixel offset by requested step size
  p.sensor.pixel_i += stepi;
  view_cam_sample(&p);
  // intersect with plane at v[1]:
  for(int k=0;k<3;k++) o[k] = path->v[1].hit.x[k] - p.v[0].hit.x[k];
  dotpo = dotproduct(path->v[1].hit.n, o);
  dotpd = dotproduct(path->v[1].hit.n, p.e[1].omega);
  t = dotpo/dotpd;
  for(int k=0;k<3;k++) rd_i[k] = p.v[0].hit.x[k] + t*p.e[1].omega[k] - path->v[1].hit.x[k];
  if(!(dotproduct(rd_i,rd_i) >= 0.0)) return 1;

  // sample new pixel offset by requested step in j direction
  p.sensor.pixel_i = path->sensor.pixel_i;
  p.sensor.pixel_j += stepj;
  view_cam_sample(&p);
  // intersect with plane at v[1]:
  for(int k=0;k<3;k++) o[k] = path->v[1].hit.x[k] - p.v[0].hit.x[k];
  dotpo = dotproduct(path->v[1].hit.n, o);
  dotpd = dotproduct(path->v[1].hit.n, p.e[1].omega);
  t = dotpo/dotpd;
  for(int k=0;k<3;k++) rd_j[k] = p.v[0].hit.x[k] + t*p.e[1].omega[k] - path->v[1].hit.x[k];
  if(!(dotproduct(rd_j,rd_j) >= 0.0)) return 1;
  return 0;
}

// computes 2x2 matrices R_k : x_1 |--> h_k
// ray differentials only make sense when actually mutating the segment from
// camera, but this also computes the transfer matrices Tp, indicating how do
// vertices change wrt to p->v[s+1].
// return dh_dx
static inline double raydifferentials_compute_rd_h(
    path_t *p,
    float *R,
    const int s,  // first vertex
    const int e)  // last vertex
{
  if(e - s <= 1)
  { // no half vec segment, no inner vertices.
    return 1.0;
  }

  const int n = e-s;
  float Li[n+1][9];
  float A[n+1][9];
  float u[n+1][9];
  float T[n+1][9];
  float H[n+1][9];
  float tmp[9], tmp2[9];

  // first pass: store Li and A

  double dh_dx = 1.0;
  if((dh_dx = _rd_inv(p, s, p->v[s].diffgeo.b, Li[0])) == 0.0f) return 0.0;
  for(int k=1;k<n;k++)
  {
    mat3_mul(p->v[s+k].diffgeo.a, Li[k-1], A[k]);
    mat3_mul(Li[k-1], p->v[s+k-1].diffgeo.c, u[k-1]);
    mat3_mul(A[k], p->v[s+k-1].diffgeo.c, tmp);
    mat3_sub(p->v[s+k].diffgeo.b, tmp, tmp2);
    // L = b - au
    if((dh_dx *= _rd_inv(p, s+k, tmp2, Li[k])) == 0.0f) return 0.0;
  }

  // loop over all inner vertex h matrices:
  if(R) for(int hi=1;hi<n;hi++)
  {
    // now finding the block x1 -> h[hi].
    // we pretend we have some input matrices h[.] = 0 except h[hi] = I.
    memset(H, 0, sizeof(float)*9*hi);
    mat3_set_identity(H[hi]);

    for(int k=hi+1;k<n;k++)
    {
      mat3_mul(A[k], H[k-1], tmp);
      for(int i=0;i<9;i++) H[k][i] = -tmp[i];
    }

    float x[9];
    mat3_mul(Li[n-1], H[n-1], x); // x[n-1]
    // TODO: can we potentially put the matrix upside down and reach this in one step always?
    for(int k=n-2;k>=1;k--)
    {
      mat3_mul(p->v[s+k].diffgeo.c, x, tmp); // x[k+1]
      for(int i=0;i<9;i++) tmp[i] = H[k][i] - tmp[i];
      mat3_mul(Li[k], tmp, x); // x[k]
    }
    // now x is the matrix that maps h[hi] -> x[1].
    // we need to invert that.
    // for volume vertices, the last row will be 0.
    // x[8] = 1.0; // XXX meaningful? analogy to R[1] = b[1] is given at least.
    // XXX FIXME: this gloriously fails for 4vtx paths in media!
    // XXX FIXME: when trying to invert the 3x3 block (2x2 seems to work, but is of course slightly meaningless)
    if(_rd_inv(p, s+hi, x, R + 9*(s+hi)) == 0.0f)
    {
      mat3_set_identity(R + 9*(s+hi));
      // fprintf(stderr, "piece of crap map h[%d] -> x[1]:\n", s+hi);
      // fprintf(stderr, "%g %g %g\n", x[0], x[1], x[2]);
      // fprintf(stderr, "%g %g %g\n", x[3], x[4], x[5]);
      // fprintf(stderr, "%g %g %g\n", x[6], x[7], x[8]);
      // path_print_manifold_matrix(p, stdout);
      // return 0.0f;
    }
    // note that for a three-vertex path the only inner matrix R[1] = b[1].
  }

  // also compute transfer matrix T_1 for measurement contribution computation
  mat3_mul(Li[n-1], p->v[e-1].diffgeo.c, T[n-1]);
  for(int k=0;k<9;k++) T[n-1][k] = -T[n-1][k];
  p->v[e-1].diffgeo.Tp[0] = T[n-1][0];
  p->v[e-1].diffgeo.Tp[1] = T[n-1][1];
  p->v[e-1].diffgeo.Tp[2] = T[n-1][3];
  p->v[e-1].diffgeo.Tp[3] = T[n-1][4];
  for(int i=n-2;i>0;i--)
  {
    mat3_mul(u[i], T[i+1], T[i]);
    for(int k=0;k<9;k++) T[i][k] = -T[i][k];
    p->v[s+i].diffgeo.Tp[0] = T[i][0];
    p->v[s+i].diffgeo.Tp[1] = T[i][1];
    p->v[s+i].diffgeo.Tp[2] = T[i][3];
    p->v[s+i].diffgeo.Tp[3] = T[i][4];
  }
  // more precise path for 3 vertex case
  if(n == 2 && R)
    for(int k=0;k<9;k++) R[9*(s+1)+k] = p->v[s+1].diffgeo.b[k];
  return fabs(dh_dx);
}

// compute the (isotropic, radius of the) pixel footprint at the given path vertex.
static inline int raydifferentials_vx(
    const path_t *path,  // path
    const int v,         // vertex index where the ray differential chain ends
    const int use_bsdf,  // flag to include bsdf roughness estimate
    float *R)            // 2x2 matrix mapping pixel to vertex area displacement in dpdu/dpdv at v[v]
{
  assert(v > 0);
  assert(v < path->length);

  // compute offsets in tangent space at vertex 1:
  float rd_i[3] = {0}, rd_j[3] = {0};
  const int err = raydifferentials_v1(path, 1.0, 1.0, rd_i, rd_j);
  if(err) return 1; // out of camera

  // express in tangent frame of v[1]
  R[0] = dotproduct(path->v[1].diffgeo.dpdu, rd_i);
  R[1] = dotproduct(path->v[1].diffgeo.dpdu, rd_j);
  R[2] = dotproduct(path->v[1].diffgeo.dpdv, rd_i);
  R[3] = dotproduct(path->v[1].diffgeo.dpdv, rd_j);

  // transform to given vertex using manifold derivatives:
  for(int i=1;i<v;i++)
  {
    // TODO: since we can't invert c, maybe it's more clever to:
    // TODO: map dx_i -> dh_{i+1} via a_{i+1}
    // TODO: map dh_{i+1} -> dx_{i+1} via b^{-1}_{i+1}
    // TODO: only that now in the last step we don't have a,b,c_{v} :(
    // TODO: work around it by going dh_{v-1} -> dx_{v} via c^{-1}_{v-1} ?
    float tmp[4], ci[4];
    // convert from dx_i to dh_i:
    const float b[4] = {
      path->v[i].diffgeo.b[0],
      path->v[i].diffgeo.b[1],
      path->v[i].diffgeo.b[3],
      path->v[i].diffgeo.b[4]};
    // FIXME: envmaps (don't have dpdu/dpdv?)!!!
    const float c[4] = {
      path->v[i].diffgeo.c[0],
      path->v[i].diffgeo.c[1],
      path->v[i].diffgeo.c[3],
      path->v[i].diffgeo.c[4]};
    mat2_mul(b, R, tmp);
    // multiply by roughness:
    if(use_bsdf)
      // TODO: check what happens in volumes (also whether we need/want 3x3 at all)
      for(int k=0;k<4;k++) tmp[k] *= path->v[i].shading.roughness;
    // convert from dh_i to dx_{i+1}
    // fprintf(stderr, "c[%d/%d] = %g %g -- %g %g \n", i, v,
        // c[0], c[1], c[2], c[3]);
    if(mat2_invert(c, ci) == 0.0f) return 1+i;
    mat2_mul(tmp, ci, R);
  }
  return 0;
}

#if 0
static inline double raydifferentials_compute_rd_h_double(
    path_t *p,
    double *R,
    const int s,  // first vertex
    const int e)  // last vertex
{
  if(e - s <= 1)
  { // no half vec segment, no inner vertices.
    return 1.0;
  }

  const int n = e-s;
  double Li[n+1][4];
  double A[n+1][4];
  double H[n+1][4];
  double tmp[4], tmp2[4];

  // new temp vars in double, too:
  double a[n+1][4];   // 3 input
  double b[n+1][4];
  double c[n+1][4];
  double u[n+1][4];   // 2 output
  double Tp[n+1][4];

  for(int k=0;k<=n;k++)
  {
    for(int i=0;i<4;i++) a[k][i] = p->v[s+k].diffgeo.a[i];
    for(int i=0;i<4;i++) b[k][i] = p->v[s+k].diffgeo.b[i];
    for(int i=0;i<4;i++) c[k][i] = p->v[s+k].diffgeo.c[i];
  }

  // first pass: store Li and A

  double dh_dx = 1.0;
  if((dh_dx = mat2d_invert(b[0], Li[0])) == 0.0f) return 1;
  for(int k=1;k<n;k++)
  {
    mat2d_mul(a[k], Li[k-1], A[k]);
    mat2d_mul(Li[k-1], c[k-1], u[k-1]);
    mat2d_mul(A[k], c[k-1], tmp);
    mat2d_sub(b[k], tmp, tmp2);
    // L = b - au
    if((dh_dx *= mat2d_invert(tmp2, Li[k])) == 0.0f) return 1;
  }
  p->cache.dh_dx = fabs(dh_dx);

  // loop over all inner vertex h matrices:
  for(int hi=1;hi<n;hi++)
  {
    // now finding the block x1 -> h[hi].
    // we pretend we have some input matrices h[.] = 0 except h[hi] = I.
    memset(H, 0, sizeof(double)*4*hi);
    mat2d_set_identity(H[hi]);

    for(int k=hi+1;k<n;k++)
    {
      mat2d_mul(A[k], H[k-1], tmp);
      for(int i=0;i<4;i++) H[k][i] = -tmp[i];
    }

    double x[4];
    mat2d_mul(Li[n-1], H[n-1], x); // x[n-1]
    for(int k=n-2;k>=1;k--)
    {
      mat2d_mul(c[k], x, tmp); // x[k+1]
      for(int i=0;i<4;i++) tmp[i] = H[k][i] - tmp[i];
      mat2d_mul(Li[k], tmp, x); // x[k]
    }
    // now x is the matrix that maps h[hi] -> x[1].
    // we need to invert that:
    if(R) if(mat2d_invert(x, R+4*(s+hi)) == 0.0f) return 1;
    // note that for a three-vertex path the only inner matrix R[1] = b[1].
  }

  // also compute transfer matrix T_1 for measurement contribution computation
  // (and wenzel's reprojection)
  mat2d_mul(Li[n-1], c[n-1], Tp[n-1]);
  for(int k=0;k<4;k++) Tp[n-1][k] = -Tp[n-1][k];
  for(int i=n-2;i>0;i--)
  {
    mat2d_mul(u[i], Tp[i+1], Tp[i]);
    for(int k=0;k<4;k++) Tp[i][k] = -Tp[i][k];
  }
#if 1 // fast path for 3 vertex case
  if(n == 2)
  {
    if(R) for(int k=0;k<4;k++) R[4*(s+1)+k] = p->v[s+1].diffgeo.b[k];
  }
#endif
  // copy back results to float matrices:
  for(int k=0;k<=n;k++)
    for(int i=0;i<4;i++) p->v[s+k].diffgeo.Tp[i] = Tp[k][i];

  return 0;
}
#endif

#endif
