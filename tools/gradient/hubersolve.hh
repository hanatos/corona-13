#pragma once
#include <cmath>
#include <Eigen/SparseCholesky>
#include <cstdint>

using namespace Eigen;

void huber_solve(
    float *primal,
    float *grad_x,
    float *grad_y,
    uint64_t wd,
    uint64_t ht,
    float alpha) // primal trust
{
  const uint64_t size = wd*ht;
  const int maxit = 7;
  for(int c=0;c<3;c++)
  {
    fprintf(stderr, "[huber] solving channel %d", c);
    VectorXd x(size);
    VectorXd w(3*size);
    VectorXd w2(3*size);
    for(uint64_t i=0;i<size;i++)
    {
      w(i) = alpha;
      w(size+i) = 1.0;
      w(2*size+i) = 1.0;
    }
    VectorXd r(3*size);
    VectorXd b(3*size);
    VectorXd h(3*size);
    SparseMatrix<double> A(3*size, size);
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(5*size);
    for(uint64_t j=0;j<ht;j++) for(uint64_t i=0;i<wd;i++)
    {
      uint64_t k = j*wd+i;
      triplets.push_back(T(k,k,w(k)));
      if(i<wd-1)
      {
        triplets.push_back(T(size+k,k, -w(k)));
        triplets.push_back(T(size+k,k+1,w(k+1)));
      }
      if(j<ht-1)
      {
        triplets.push_back(T(2*size+k,k,  -w(k)));
        triplets.push_back(T(2*size+k,k+wd,w(k+wd)));
      }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    for(int it=0;it<maxit;it++)
    {
      fprintf(stderr, ".");
      fflush(stderr);
      // fill b
      for(int i=0;i<size;i++)
      {
        b(i)        = w(i)        * primal[3*i+c];
        b(i+size)   = w(i+size)   * grad_x[3*i+c];
        b(i+2*size) = w(i+2*size) * grad_y[3*i+c];
      }
      // solve Qx = h
      // (A' W^2 A) x = A' W^2 b
      h = A.transpose() * b;
      SimplicialLDLT<SparseMatrix<double> > solver;
      SparseMatrix<double> Q = SparseMatrix<double>(A.transpose()) * A;
      solver.compute(Q);
      if(solver.info()!=Success)
      {
        fprintf(stderr, "aeaoaarrgh\n");
        // decomposition failed
        return;
      }
      x = solver.solve(h);
      if(solver.info()!=Success)
      {
        fprintf(stderr, "sthrchsrch\n");
        // solving failed
        return;
      }

      // residual:
      r = A * x - b;

      // eps = max|b|/100
      const double eps = 0.01;

      // select weights based on huber norm:
      // w_ii = 1/(1 + (r_i/eps)^2)^1/4
      for(int i=0;i<3*size;i++)
      {
        w2(i) = 1.0/std::pow(1.0 + (r(i)/eps)*(r(i)/eps), 0.25);
        if(i < size) w2(i) *= alpha;
      }

      for (int k=0; k<A.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
          A.coeffRef(it.row(), it.col()) *= w2(it.col())/w(it.col());

      for(int i=0;i<3*size;i++)
      {
        w(i) = 1.0/std::pow(1.0 + (r(i)/eps)*(r(i)/eps), 0.25);
        if(i < size) w(i) *= alpha;
      }
    }
    for(uint64_t k=0;k<size;k++)
      primal[3*k+c] = x(k);
    fprintf(stderr, "\n");
  }
}
