#ifndef _FINITE_DIFF_H_
#define _FINITE_DIFF_H_

#include "eigen/Eigen/Core"
#include <functional>
#include <iostream>
#include <cmath>

static const double eps = 1e-3;

using Eigen::VectorXd;
using Eigen::MatrixXd;

double finite_diff_scalar2scalar(std::function<double(double)> f, double x);

VectorXd finite_diff_vec2scalar(std::function<double(VectorXd)> f, VectorXd x);

MatrixXd finite_diff_vec2vec(std::function<VectorXd(VectorXd)> f, VectorXd x, int out_size);

MatrixXd finite_diff2_vec2scalar(std::function<double(VectorXd, VectorXd)> f, VectorXd x1, VectorXd x2);

#endif
