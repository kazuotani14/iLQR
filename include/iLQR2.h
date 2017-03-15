#pragma once
#include "standardIncludes.h"

#define epsilon 1e-3

class iLQR
{
  Eigen::VectorXd* Alpha;
  const double tolFun;
  const double tolGrad;
  const int maxIter;
  const double lambda;
  const double dlambda;
  const double lambdaFactor;
  const double lambdaMax;
  const double lambdaMin;
  const double zMin;

  const int n;          // dimension of state vector
  const int m;          // dimension of control vector
  const int T;          // number of state transitions (#timesteps-1)

  std::vector<Eigen::VectorXd> xs;
  std::vector<Eigen::VectorXd> us;
  std::vector<Eigen::VectorXd> ls;
  std::vector<Eigen::MatrixXd> Ls;
  Eigen::VectorXd xd; // target state

  Eigen::Matrix2d control_limits;
  Eigen::VectorXd u0; //iniital control sequence
  // trace, for debugging? struct?

  // Helper functions will modify member variables instead of returning values,
  //  because most of them need to return multiple values.
  void forward_pass(Eigen::VectorXd &x0, Eigen::Matrix2d &u, Eigen::MatrixXd &L,
                    Eigen::Matrix2d &x, Eigen::Matrix2d &du, Eigen::VectorXd &alpha);
  // void back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,regType,lims,u);
  // void boxQP(H,g,lower,upper,x0,options);


public:
  iLQR(): tolFun(pow(10,-5)), tolGrad(pow(10,-5)), maxIter(30), lambda(1),
                dlambda(1), lambdaFactor(1.6), lambdaMax(pow(10,10)), lambdaMin(pow(10,-7)),
                zMin(0), n(6), m(2), T(50)
  {
    // TODO find way to do this without pointer
    Alpha = new Eigen::VectorXd(11);
    *Alpha << 1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    // TODO initialize control inputs here
  }

  ~iLQR(){
    delete Alpha;
  }

  void generate_trajectory(const Eigen::VectorXd &x_0, const Eigen::VectorXd &x_d, int trajectoryLength, int nControlInputs);

};
