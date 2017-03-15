#ifndef _ILQR_H_
#define _ILQR_H_

#include "standardIncludes.h"

#define epsilon 1e-3

class iLQR
{
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
  Eigen::VectorXd Alpha;

  // Tracking progress
  std::vector<Eigen::Matrix2d> xs; // states from last trajectory
  std::vector<Eigen::Matrix2d> us; // controls from last trajectory
  std::vector<Eigen::VectorXd> ls;
  std::vector<Eigen::MatrixXd> Ls;

  // Stuff relevant at each cycle
  Eigen::Matrix2d control_limits;
  Eigen::Matrix2d u0; //initial control sequence
  Eigen::VectorXd x_current;
  Eigen::VectorXd u_current;

  // Helper functions will modify passed-in variables instead of returning values,
  //  because most of them need to return multiple values.
  // void init_traj(Eigen::VectorXd &x_0, Eigen::Matrix2d &u0);
  // void forward_pass(Eigen::VectorXd &x0, Eigen::Matrix2d &u, Eigen::MatrixXd &L,
  //                   Eigen::Matrix2d &x, Eigen::Matrix2d &du, Eigen::VectorXd &alpha);
  // int back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,regType,lims,u); //TODO
  // void boxQP(H,g,lower,upper,x0);
	//function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = car_dyn_cst(x,u,full_DDP)

public:
  Eigen::VectorXd x_d; // target state

  iLQR(): tolFun(pow(10,-5)), tolGrad(pow(10,-5)), maxIter(30), lambda(1),
                dlambda(1), lambdaFactor(1.6), lambdaMax(pow(10,10)), lambdaMin(pow(10,-7)),
                zMin(0), n(6), m(2), T(50)
  {
    // TODO figure out how to do this
    // Alpha << 1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    // TODO read control inputs from somewhere
    // init_traj(x_0, u0);
  }

  // void generate_trajectory(const Eigen::VectorXd &x_0, const Eigen::VectorXd &x_d, int trajectoryLength, int nControlInputs);

  double get_dynamics_and_cost(const Eigen::VectorXd &x, const Eigen::VectorXd u, Eigen::VectorXd &dx);

  // Pure virtual functions, to be overridden in LocoCar
  virtual Eigen::VectorXd dynamics(const Eigen::VectorXd &x, const Eigen::Vector2d &u) = 0; //dynamics
  virtual double cost(const Eigen::VectorXd &x, const Eigen::VectorXd &u) = 0;
  virtual double final_cost(const Eigen::VectorXd &x) = 0;
};

#endif
