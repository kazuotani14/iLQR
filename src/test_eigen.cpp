#include <iostream>
#include <eigen/Eigen/Dense>
#include <math.h>

using namespace Eigen;

template<typename T>
double first_value(T &mat){
  return mat(0);
}

VectorXd first_val_eigen(VectorXd &mat){
  return mat;
}

double pi = M_PI;


int main()
{

  // double alpha_F, alpha_R;

  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);

  VectorXd v;
  v << 1, 2, 3, 4, 5;
  //v.head(3), v.tail(2)
  //std::abs(1.3)

  std::cout << first_val_eigen(v) << std::endl;
}
