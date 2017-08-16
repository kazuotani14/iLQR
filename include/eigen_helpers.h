#ifndef _EIGEN_HELPERS_H_
#define _EIGEN_HELPERS_H_

#include "common.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;


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


//TODO make this faster, with no resizes, or use std
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

//TODO there has to be a better way to do this
// Equivalent of `mat = mat(bool_vec, bool_vec)` in matlab
inline MatrixXd extract_bool_rowsandcols(const MatrixXd& mat, const VectorXd& bool_vec)
{
  int n_dims = bool_vec.sum();
  MatrixXd small_mat(n_dims, n_dims);

  int row_idx = 0;
  for(int i=0; i<mat.rows(); i++)
  {
    int col_idx = 0;
    if(bool_vec(i) != 0)
    {
      for(int j=0; j<mat.cols(); j++)
      {
        if(bool_vec(j) != 0) small_mat(row_idx, col_idx++) = mat(i,j);
      }
      row_idx++;
    }
  }
  return small_mat;
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


inline VectorXd elem_square(const VectorXd &vec)
{
  return vec.array().square().matrix();
}

inline VectorXd elem_sqrt(const VectorXd &vec)
{
  return vec.array().sqrt().matrix();
}

// Differentiable "soft" absolute value function
inline VectorXd sabs(const VectorXd &vec, const VectorXd &p)
{
  VectorXd sum = elem_sqrt(elem_square(vec)+elem_square(p));
  return sum - p;
}

#endif
