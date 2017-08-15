#ifndef _BOXQP_H_
#define _BOXQP_H_

#include "eigen_helpers.h"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// #define VERBOSE

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Optimization parameters
static const int qp_maxIter        = 100;      // maximum number of iterations
static const double minGrad        = 1e-8;      // minimum norm of non-fixed gradient
static const double minRelImprove  = 1e-8;      // minimum relative improvement
static const double stepDec        = 0.6;       // factor for decreasing stepsize
static const double minStep        = 1e-22;     // minimal stepsize for linesearch
static const double Armijo         = 0.1;  // Armijo parameter (fraction of linear improvement required in line search)

struct lineSearchResult
{
  lineSearchResult(int n_dims): x_opt(n_dims), v_opt(n_dims) {}

  bool failed = false;
  int n_steps = 0;
  VectorXd x_opt;
  double v_opt;
};

struct boxQPResult
{
  boxQPResult(int n_dims):
    x_opt(n_dims), v_free(n_dims), H_free(n_dims, n_dims) {}

  int result = 0;
  VectorXd x_opt;
  VectorXd v_free;
  MatrixXd H_free;
};

VectorXd clamp_to_limits(const VectorXd &x, const VectorXd& lower, const VectorXd& upper);

double quadCost(const MatrixXd& Q, const VectorXd& c, const VectorXd& x);

lineSearchResult quadclamp_line_search(const VectorXd& x0, const VectorXd& search_dir,
                     const MatrixXd& Q, const VectorXd& c,
                     const VectorXd& lower, const VectorXd& upper);


boxQPResult boxQP(const MatrixXd &Q, const VectorXd &c, const VectorXd &x0,
                  const VectorXd& lower, const VectorXd& upper);

inline bool approx_eq(double a, double b)
{
  return (std::abs(a-b) < 1e-4);
}

#endif
