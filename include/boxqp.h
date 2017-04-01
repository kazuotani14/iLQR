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

/*
	Minimize 0.5*x'*H*x + x'*g  s.t. lower<=x<=upper

	 inputs:
	    H            - positive definite matrix  (m * m)
	    g            - bias vector         			 (m)
	    x0           - initial state             (m)

			lower        - lower bounds              (m)
	    upper        - upper bounds              (m)

	 outputs:
	 		result       - result type (roughly, higher is better, see below)
	    x            - solution = k_i             (m)
	    res.H_free        - subspace cholesky factor=R (n_free * n_free)
	    res.v_free         - set of free dimensions     (m)
											- vector of 0 or 1
*/

// TODO rename (H,g) to (Q,c)

// Optimization parameters
static const int qp_maxIter           = 100;      // maximum number of iterations
static const double minGrad        = 1e-8;      // minimum norm of non-fixed gradient
static const double minRelImprove  = 1e-8;      // minimum relative improvement
static const double stepDec        = 0.6;     	// factor for decreasing stepsize
static const double minStep        = 1e-22;     // minimal stepsize for linesearch
static const double Armijo         = 0.1;   		// Armijo parameter (fraction of linear improvement required in line search)

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

double quadCost(const MatrixXd& H, const VectorXd& g, const VectorXd& x);

VectorXd clamp_to_limits(const VectorXd &x, const VectorXd& lower, const VectorXd& upper);

lineSearchResult quadclamp_line_search(const VectorXd x0, const VectorXd search_dir,
														 const MatrixXd H, const VectorXd g,
														 const VectorXd lower, const VectorXd upper);


boxQPResult boxQP(const MatrixXd &H, const VectorXd &g, const VectorXd &x0,
                  const VectorXd& lower, const VectorXd& upper);
				  
#endif
