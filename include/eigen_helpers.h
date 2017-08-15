#ifndef _EIGEN_HELPERS_H_
#define _EIGEN_HELPERS_H_

#include "common.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

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

// Extracts elements of vec for which indices is non-zero
inline VectorXd subvec_w_ind(const VectorXd& vec, const VectorXd& indices)
{
	assert(vec.size() == indices.size());
	std::vector<double> vals;
	for(int i=0; i<vec.size(); i++){
		if(indices(i)>0){
			vals.push_back(vec(i));
		}
	}
	Eigen::Map<VectorXd> subvec(vals.data(), vals.size());
	return subvec;
}

//TODO make this faster, with no resizes
inline MatrixXd rows_w_ind(MatrixXd &mat, VectorXd &indices)
{
  MatrixXd submat;
  if(mat.rows() != indices.size()){
    std::cout << "mat.rows != indices.size\n";
    return submat;
  }
	for(int i=0; i<indices.size(); i++){
		if(indices(i)>0){
      submat.conservativeResizeLike(MatrixXd(submat.rows()+1, mat.cols()));
			submat.row(submat.rows()-1) = mat.row(i);
		}
	}
	return submat;
}

const Eigen::IOFormat CleanFmt(3, 0, ", ", "\n", "", "");

inline void print_eigen(const std::string name, const Eigen::Ref<const Eigen::MatrixXd>& mat)
{
  if(mat.cols() == 1)
  {
    std::cout << name << ": " << mat.transpose().format(CleanFmt) << ";" << std::endl;
  }
  else{
    std::cout << name << ":\n" << mat.format(CleanFmt) << std::endl;
  }
}

#endif
