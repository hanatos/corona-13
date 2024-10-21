#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// http://www3.cplusplus.com/forum/general/222617/
#define real float

void madd(
    real *const restrict u,
    const float m,
    const real *const restrict v,
    const uint64_t size)
{
#pragma omp parallel for simd schedule(static)
  for(int i=0;i<size;i++)
    u[i] += m*v[i];
}

void addm(
    const real *const restrict u,
    const float m,
    real *const restrict v,
    const uint64_t size)
{
#pragma omp parallel for simd schedule(static)
  for(int i=0;i<size;i++)
    v[i] = m*v[i] + u[i];
}

real cg_dot(const real *const restrict u, const real *const restrict v, const uint64_t size)
{
  real res = 0.0;
#pragma omp parallel for schedule(static) reduction(+:res) default(shared)
  for(int i=0;i<size;i++)
    res += u[i]*v[i];
  return res;
}

// solve Ax = b where x is the output image, b : [primal, grad_x, grad_y] and A contains the 2D poisson process problem
// actually solve Qx = h with weights W such that
// (A' W^2 A) x = A' W^2 b
// TODO: solve on colours directly, the whole matrix/cg dance is exactly the same
void cg_solve(
    float *primal,
    float *grad_x,
    float *grad_y,
    uint64_t wd,
    uint64_t ht,
    float in_alpha) // primal trust
{
  const real cg_eps = 1e-5;

  const uint64_t size = wd*ht;
  // primal output, pixel sized
  real *x = malloc(sizeof(real)*size);
  // temporary result, pixel sized
  real *p = malloc(sizeof(real)*size);
  // cached (Q p) pixels
  real *Qp = malloc(sizeof(real)*size);
  // residual h - Q p in pixels
  real *r = malloc(sizeof(real)*size);
  // error b - A x including gradients
  real *e = malloc(sizeof(real)*size*3);
  real *b = malloc(sizeof(real)*size*3);
  // weights needed to construct matrix wA
  real *w = malloc(sizeof(real)*size*3);

  for(int c=0;c<3;c++)
  {
    // init weights:
#pragma omp parallel for schedule(static)
    for(uint64_t k=0;k<size;k++)
    {
      w[k] = 1.0;
      w[k+size] = 1.0;
      w[k+2*size] = 1.0;
    }
    memset(x, 0, sizeof(real)*size);
#pragma omp parallel for schedule(static)
    for(uint64_t k=0;k<size;k++)
    {
      b[k]        = in_alpha * primal[3*k+c];
      b[k+size]   = grad_x[3*k+c];
      b[k+2*size] = grad_y[3*k+c];
    }

    // huber norm iteration/weighting:
    for(int outerit=0;outerit<20;outerit++)
    {
      // memset(x, 0, sizeof(real)*size);
      // init error e = b - Ax
#pragma omp parallel for schedule(static) collapse(2)
      for(uint64_t j=0;j<ht;j++) for(uint64_t i=0;i<wd;i++)
      {
        const uint64_t k = j*wd+i;
        e[k]        = b[k] - in_alpha * x[k];
        e[k+  size] = i < wd-1 ? b[k+  size] - (x[k+1]  - x[k]) : 0.0;
        e[k+2*size] = j < ht-1 ? b[k+2*size] - (x[k+wd] - x[k]) : 0.0;
      }

      // update weights according to residual
      if(outerit)
      {
        const real eps = 0.0001;
        // real sum = 0.0;
#pragma omp parallel for schedule(static) //reduction(+:sum)
        for(uint64_t k=0;k<3*size;k++)
        {
          // w[k] = 1.0/pow(1.0 + (e[k]/eps)*(e[k]/eps), 0.25);
          w[k] = 1.0/sqrt(eps + e[k]*e[k]);
          // w[k] = 1.0/pow(eps + e[k]*e[k], 0.7);//0.25);
          // w[k] = 1.0;
          // w[k] = 1.0/(eps + e[k]*e[k]);
          // if(k < size) w[k] = 1.0; // L2 for primal XXX need to make sure weights are about in same range
          // w[k] = 1.0/(eps + fabs(e[k]));
          // w[k]*=w[k];
          // sum += w[k];
        }
// #pragma omp parallel for schedule(static)
        // for(uint64_t k=0;k<3*size;k++) w[k] *= size/sum;
      }
      pfm_t ptmp = {.pixel=w,.wd=wd,.ht=ht};
      pfm_write(&ptmp, "weights.pfm");

      // init residual r = h - Qx = Atw2e
#pragma omp parallel for schedule(static) collapse(2)
      for(uint64_t j=0;j<ht;j++) for(uint64_t i=0;i<wd;i++)
      {
        const uint64_t k = j*wd+i;
        r[k] = e[k] * w[k] * in_alpha;
        if(i)        r[k] += e[k+  size- 1] * w[k+  size- 1];
        if(j)        r[k] += e[k+2*size-wd] * w[k+2*size-wd];
        if(i < wd-1) r[k] -= e[k+  size]    * w[k+  size];
        if(j < ht-1) r[k] -= e[k+2*size]    * w[k+2*size];
      }

#pragma omp parallel for schedule(static)
      for(uint64_t k=0;k<size;k++)
        p[k] = r[k];

      real rr = cg_dot(r, r, size);

      // cg iteration:
      for(int innerit=0;innerit<60;innerit++)
      {
        // compute Qp
#pragma omp parallel for schedule(static) collapse(2)
        for(uint64_t j=0;j<ht;j++) for(uint64_t i=0;i<wd;i++)
        {
          const uint64_t k = j*wd+i;
          Qp[k] = p[k] * w[k] * in_alpha * in_alpha;
          if(i) Qp[k] += (p[k] - p[k- 1]) * w[k+  size- 1];
          if(j) Qp[k] += (p[k] - p[k-wd]) * w[k+2*size-wd];
          if(i < wd-1) Qp[k] += (p[k] - p[k+ 1]) * w[k+  size];
          if(j < ht-1) Qp[k] += (p[k] - p[k+wd]) * w[k+2*size];
        }

        const real alpha = rr / fmax(cg_dot(p, Qp, size), cg_eps);
        madd(x, alpha, p, size);
        madd(r, -alpha, Qp, size);

        real new_rr = cg_dot(r, r, size);
        // fprintf(stderr, "[it %d] alpha= %g res %g -> %g\n", innerit, alpha, rr, new_rr);
        if(new_rr < cg_eps) break;

        const real beta = new_rr / fmax(rr, cg_eps);
        addm(r, beta, p, size);
        rr = new_rr;
      }
    }

    // overwrite primal with x
#pragma omp parallel for schedule(static)
    for(uint64_t k=0;k<size;k++) primal[3*k+c] = x[k];
  }

  free(x);
  free(p);
  free(Qp);
  free(e);
  free(r);
  free(w);
}
