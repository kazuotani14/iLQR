#ifndef _ILQR_H_
#define _ILQR_H_

#include "standardIncludes.h"

class iLQR
{
  const double tolFun;
  const double tolGrad;
  const int maxIter;
  double lambda;
  double dlambda;
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
  double cost_s;

  // Pure virtual functions, to be overridden in LocoCar
  virtual VecXd dynamics(const VecXd &x, const VecXd &u) = 0; //dynamics
  virtual VecXd integrate_dynamics(const VecXd &x, const VecXd u) = 0;
  virtual double cost(const VecXd &x, const VecXd &u) = 0;
  virtual double final_cost(const VecXd &x) = 0;

  // Helper functions will modify passed-in variables instead of returning values,
  //  because most of them need to return multiple values.
  void forward_pass(const VecXd &x0, const VecOfVecXd &u,
  									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost);

  void forward_pass(const VecXd &x0, const VecOfVecXd &u,
  									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
  									const VecOfVecXd &x, const VecOfMatXd &L);

  int backward_pass(const VecOfVecXd &cx, const VecOfVecXd &cu, const VecOfMatXd &cxx, const VecOfMatXd &cxu,
  									const VecOfMatXd &cuu, const VecOfMatXd &fx, const VecOfMatXd &fu, const VecOfVecXd &u,
  									VecOfVecXd &Vx, VecOfMatXd &Vxx, VecOfVecXd &k, VecOfMatXd &K, Vec2d &dV);

  int boxQP(MatXd &H, VecXd &g, VecXd &x0,  VecXd &x, MatXd& Hfree, VecXd &free);
  VecXd clamp_to_limits(VecXd &u);

  MatXd rows_w_ind(MatXd &mat, VecXd &rows);
  VecXd subvec_w_ind(VecXd &vec, VecXd &indices);
  VecOfVecXd adjust_u(VecOfVecXd &u, VecOfVecXd &l, double alpha);

  // Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
  void compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd &fx,
  												 VecOfMatXd &fu, VecOfVecXd &cx, VecOfVecXd &cu,
  												 VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu);
  void get_dynamics_derivatives(const VecXd &x, const VecXd &u, MatXd &fx, MatXd &fu);
  void get_cost_derivatives(const VecXd &x, const VecXd &u, VecXd &cx, VecXd &cu);
  void get_cost_2nd_derivatives(const VecXd &x, const VecXd &u, MatXd &cxx, MatXd &cxu, MatXd &cuu);

  void demoQP();

public:
  VecXd x_d; // target state
  VecOfVecXd u0; //initial control sequence
  const double timeDelta; //dt for euler integration

  double init_traj(VecXd &x_0, VecOfVecXd &u0);

  iLQR(): tolFun(pow(10,-5)), tolGrad(pow(10,-5)), maxIter(30), lambda(1),
                dlambda(1), lambdaFactor(1.6), lambdaMax(pow(10,10)), lambdaMin(pow(10,-7)),
                zMin(0), n(6), m(2), T(50), timeDelta(0.05)
  {
 		Alpha.resize(11);
    Alpha << 1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    // TODO move this somewhere else
    control_limits.resize(2,2);
    control_limits << -1, 4,        //min/max for throttle
                      -0.76, 0.68;  //min/max for steering
  }

  // Generate and solve new trajectory optimization problem
  void generate_trajectory(const VecXd &x_0, const VecXd &x_d, int trajectoryLength);

};

#endif
