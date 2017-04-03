#ifndef _ILQR_H_
#define _ILQR_H_

#include "gtest/gtest_prod.h"
#include "standardIncludes.h"

#include "model.h"
#include "boxqp.h"
#include <memory>

static const double tolFun = 1e-6;
static const double tolGrad = 1e-6;
static const int maxIter = 30;
static double lambda = 1;
static double dlambda = 1;
static const double lambdaFactor = 1.6;
static const double lambdaMax = 1e11;
static const double lambdaMin = 1e-8;
static const double zMin = 0;
static std::vector<double> alpha_vec = {1.0000, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010};
static Eigen::Map<VectorXd> Alpha(alpha_vec.data(), alpha_vec.size());

class iLQR
{
public:
  iLQR(Model* p_dyn, double timeDelta, int timesteps): dt(timeDelta), T(timesteps)
  {
    model.reset(p_dyn);
  }
  iLQR() = default;

  VectorXd x_d; // target state
  std::shared_ptr<Model> model;

  double init_traj(VectorXd &x_0, VecOfVecXd &u0);
  void generate_trajectory();

  // Tester functions
  void output_to_csv();

private:
  FRIEND_TEST(ILQRSetup, dDynamicsTest);
  FRIEND_TEST(ILQRSetup, dCostTest);
  FRIEND_TEST(ILQRSetup, ddCostTest);
  FRIEND_TEST(ILQRSetup, ForwardPassTest);

  int T;  // number of state transitions
  double dt;

  // Tracking progress
  VecOfVecXd xs; // states from last trajectory
  VecOfVecXd us; // controls from last trajectory
  VecOfVecXd ls;
  VecOfMatXd Ls;

  // Stuff relevant at each cycle
  VectorXd x_current;
  VectorXd u_current;
  double cost_s;

  MatrixXd du; 	//2*T double
  VecOfMatXd fx; //nxnx(T+1)
  VecOfMatXd fu; //nx2x(T+1)
  VecOfVecXd cx; //nx(T+1)
  VecOfVecXd cu; //2x(T+1)
  VecOfMatXd cxx; //nxnx(T+1)
  VecOfMatXd cxu; //nx2x(T+1)
  VecOfMatXd cuu; //2x2x(T+1)

  Vector2d dV; //2x1
	VecOfVecXd Vx; //nx(T+1)
	VecOfMatXd Vxx; //nxnx(T+1)
	VecOfMatXd L; //2xnxT
	VecOfVecXd l; //2xT

  double forward_pass(const VectorXd& x0, const VecOfVecXd& u);
  double forward_pass(const VectorXd& x0, const VecOfVecXd& u,
  									const VecOfVecXd& x, const VecOfMatXd& L);
  int backward_pass();
  VecOfVecXd adjust_u(VecOfVecXd &u, VecOfVecXd &l, double alpha);

  // Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
  void compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
};

#endif
