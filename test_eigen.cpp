#include <iostream>
#include <eigen/Eigen/Dense>
#include <math.h>
#include <vector>
#include "include/standardIncludes.h"


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

  std::cout << "K_i: " << K_i.rows() << ' ' << K_i.cols() << '\n';
  std::cout << "Quu: " << Quu.rows() << ' ' << Quu.cols() << '\n';

  MatXd test1 =  K_i.transpose()*Quu*K_i;
  MatXd test2 =  K_i.transpose()*Qux;
  // VecXd test2 =  K_i.transpose()*Qux;

  // VecXd test_push;
  // std::cout << test_push.size() << "\n\n";
  // push_back(test_push,5);
  // std::cout << test_push << "\n\n";
  // push_back(test_push,2);
  // std::cout << test_push << std::endl;


}
