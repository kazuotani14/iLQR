#ifndef _EIGEN_HELPERS_H_
#define _EIGEN_HELPERS_H_

#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Extracts elements of vec for which indices is non-zero
VectorXd subvec_w_ind(const VectorXd& vec, const VectorXd& indices);

MatrixXd rows_w_ind(MatrixXd &mat, VectorXd &indices);

const Eigen::IOFormat CleanFmt(3, 0, ", ", "\n", "", "");

void print_eigen(const std::string name, const Eigen::Ref<const Eigen::MatrixXd>& mat);

inline VectorXd elem_square(const VectorXd &vec)
{
  return vec.array().square().matrix();
}

inline VectorXd elem_sqrt(const VectorXd &vec)
{
  return vec.array().sqrt().matrix();
}

inline VectorXd sabs(const VectorXd &vec, const VectorXd &p)
{
  //Differentiable "soft" absolute value function
  VectorXd sum = elem_sqrt(elem_square(vec)+elem_square(p));
  return sum - p;
}

#endif
