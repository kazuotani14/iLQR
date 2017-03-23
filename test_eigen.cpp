#include <iostream>
#include <eigen/Eigen/Dense>
#include <math.h>
#include <vector>
#include "include/standardIncludes.h"

MatXd rows_w_ind(MatXd &mat, VecXd &indices)
{
  MatXd submat;
  if(mat.rows() != indices.size()){
    std::cout << "mat.rows != indices.size\n";
    return submat;
  }
	for(int i=0; i<indices.size(); i++){
		if(indices(i)>0){
      std::cout << i << '\n';
      submat.conservativeResizeLike(MatXd(submat.rows()+1, mat.cols()));
			submat.row(submat.rows()-1) = mat.row(i);
		}
	}
	return submat;
}

double x_squared(double x){
	double y = x*x;
  return y;
}

double finite_differences(double x0, double(*foo)(double))
{
  double epsilon = 0.00001;
  double dx = (foo(x0+epsilon)-foo(x0-epsilon)) / (2*epsilon);
  return dx;
}

VecXd vec_squared(const VecXd &x)
{
  VecXd squared = x.array().square().matrix();
  return squared;
}

MatXd finite_differences(VecXd &x,
                                    VecXd (*foo)(const VecXd&))
{
  double eps = 0.0001;
  int n = x.size();
  MatXd J(n,n);

  for (int i=0; i<n; i++)
  {
    VecXd plus = x;
    VecXd minus = x;
    plus(i) += eps;
    minus(i) -= eps;
    J.col(i) = (foo(plus)-foo(minus)) / (2*eps);
  }
  return J;
}

int main()
{

  for (int i=0; i<5; i++){
    std::cout << Vec2d::Random() << "\n\n";
  }


  Eigen::Matrix2d y;
  y << 1, 2, 3, 4;

  VecXd x(4);
  x << 1, 2, 3, 4;

  Vec2d ind;
  ind << 0, 1;

  x = x.array() / 2;

  //Testing function pointers
  VecXd (*func)(const VecXd&);
  func = &vec_squared;

  VecOfVecXd vector(3);

  int m = 2;
  int n = 6;

  VecXd Qx(n);
	VecXd Qu(m);
	MatXd Qxx(n,n);
	MatXd Qux(m,n);
	MatXd Quu(m,m);
	VecXd k_i(m);
	MatXd K_i(m,n);


  MatXd Lfree(m,n);
  MatXd R(2,2);
  R <<  1.0012, -0.0069,
         0 ,   1.5170;
  MatXd Qux_reg(m,n);
  Qux_reg << -0.0124,    0.0101,    0.0113 ,   0.0100,   -0.0196,    0.0065 ,
              -0.5332,    0.1886,   -0.7361,    1.7398,   -2.3677,   -0.4859;
  VecXd free_v(2); free_v << 1,1;

  // std::cout << rows_w_ind(Qux_reg, free_v) << "\n\n";
  // std::cout << "-------\n";
  //
  // Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, free_v));

  // std::cout << Lfree << '\n';
  // VecXd test2 =  K_i.transpose()*Qux;

  // VecXd test_push;
  // std::cout << test_push.size() << "\n\n";
  // push_back(test_push,5);
  // std::cout << test_push << "\n\n";
  // push_back(test_push,2);
  // std::cout << test_push << std::endl;


}
