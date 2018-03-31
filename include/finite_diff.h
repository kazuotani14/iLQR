#ifndef _FINITE_DIFF_H_
#define _FINITE_DIFF_H_

#include "eigen/Eigen/Core"
#include <iostream>
#include <cmath>

static const double eps = 1e-3;

using Eigen::VectorXd;
using Eigen::MatrixXd;

inline double finite_diff_scalar2scalar(std::function<double(double)> f, double x) {
  double plus = x+eps;
  double minus = x-eps;
  return (f(plus)-f(minus)) / (2*eps);
}

inline VectorXd finite_diff_vec2scalar(std::function<double(VectorXd)> f, VectorXd x) {
  int n_dims = x.size();
  VectorXd plus(n_dims), minus(n_dims), dx(n_dims);

  for (int i=0; i<n_dims; i++) {
    plus = minus = x;
    plus(i) += eps;
    minus(i) -= eps;
    dx(i) = (f(plus)-f(minus)) / (2*eps);
  }
  return dx;
}

inline MatrixXd finite_diff_vec2vec(std::function<VectorXd(VectorXd)> f, VectorXd x, int out_size) {
  int n_dims = x.size();
  VectorXd plus(n_dims), minus(n_dims);
  MatrixXd dx(out_size, n_dims);

  for (int i=0; i<n_dims; i++) {
    plus = minus = x;
    plus(i) += eps;
    minus(i) -= eps;
    dx.col(i) = (f(plus)-f(minus)) / (2*eps);
  }
  return dx;
}

inline MatrixXd finite_diff2_vec2scalar(std::function<double(VectorXd, VectorXd)> f, VectorXd x1, VectorXd x2) {
  MatrixXd dx1x2(x1.size(), x2.size());
  VectorXd p1, p2, m1, m2;

  for (int i=0; i<x1.size(); i++){
    for (int j=0; j<x2.size(); j++){
      p1 = m1 = x1;
      p2 = m2 = x2;
      p1(i) += eps;
      m1(i) -= eps;
      p2(j) += eps;
      m2(j) -= eps;
      dx1x2(i,j) = dx1x2(j,i) = (f(p1, p2) - f(m1, p2) - f(p1, m2) + f(m1, m2)) / (4*pow(eps,2));
    }
  }
  return dx1x2;
}


#endif
