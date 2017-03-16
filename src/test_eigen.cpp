#include <iostream>
#include <eigen/Eigen/Dense>
#include <math.h>
#include <vector>


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

double pi = M_PI;


int main()
{

  VecXd x(4);
  x << 1, 2, 3, 4;

  //Testing function pointers
  VecXd (*func)(const VecXd&);
  func = &vec_squared;

  VecOfVecXd vector(3);


  std::cout << finite_differences(x, func) << std::endl;
}
