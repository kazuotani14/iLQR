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
  int T;          // number of state transitions (#timesteps-1)
  VecXd Alpha;

  // Tracking progress
  VecOfVecXd xs; // states from last trajectory
  VecOfVecXd us; // controls from last trajectory
  VecOfVecXd ls;
  VecOfMatXd Ls;

  // Stuff relevant at each cycle
  MatXd control_limits;
  VecXd x_current;
  VecXd u_current;

  // Helper functions will modify passed-in variables instead of returning values,
  //  because most of them need to return multiple values.
  void forward_pass(VecXd &x0, VecOfVecXd &u,
  									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost);

  void forward_pass(VecXd &x0, VecOfVecXd &u,
  									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
  									VecOfVecXd &x, VecOfMatXd &L,
  									VecOfVecXd &du, double &alpha);

  // int back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,regType,lims,u); //TODO
  // void boxQP(H,g,lower,upper,x0);
	//function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = car_dyn_cst(x,u,full_DDP)

public:
  VecXd x_d; // target state
  VecOfVecXd u0; //initial control sequence
  const double timeDelta; //dt

  double init_traj(VecXd &x_0, VecOfVecXd &u0);

  iLQR(): tolFun(pow(10,-5)), tolGrad(pow(10,-5)), maxIter(30), lambda(1),
                dlambda(1), lambdaFactor(1.6), lambdaMax(pow(10,10)), lambdaMin(pow(10,-7)),
                zMin(0), n(6), m(2), T(50), timeDelta(0.05)
  {
    // initialize Alpha, control_limits
    //
    // TODO figure out how to do this
    // Alpha << 1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;


  }

  // Pure virtual functions, to be overridden in LocoCar
  virtual VecXd dynamics(const VecXd &x, const VecXd &u) = 0; //dynamics
  virtual double cost(const VecXd &x, const VecXd &u) = 0;
  virtual double final_cost(const VecXd &x) = 0;


  // Trajectory generation and related functions
  void generate_trajectory(const VecXd &x_0, const VecXd &x_d, int trajectoryLength);

  double get_dynamics_and_cost(const VecXd &x, const VecXd u, VecXd &dx);

  void compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd &fx,
  												 VecOfMatXd &fu, VecOfVecXd &cx, VecOfVecXd &cu,
  												 VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu);
  void get_dynamics_derivatives(const VecXd &x, const VecXd &u, MatXd &fx, MatXd &fu);
  void get_cost_derivatives(const VecXd &x, const VecXd &u, VecXd &cx, VecXd &cu);
  void get_cost_2nd_derivatives(const VecXd &x, const VecXd &u, MatXd &cxx, MatXd &cxu, MatXd &cuu);

};

#endif
