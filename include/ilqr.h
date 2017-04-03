#ifndef _ILQR_H_
#define _ILQR_H_

#include "common.h"
#include "model.h"
#include "boxqp.h"
#include <memory>
#include <numeric>

#include "gtest/gtest_prod.h"

static const int maxIter = 10;
static const double tolFun = 1e-6;
static const double tolGrad = 1e-6;
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

  std::shared_ptr<Model> model;
  VectorXd x_d; // target state TODO use this

  double init_traj(VectorXd &x_0, VecOfVecXd &u0);
  void generate_trajectory();

  void output_to_csv(const std::string filename);

private:
  int T;  // number of state transitions
  double dt;

  VectorXd x0;

  VecOfVecXd xs; // states from last trajectory
  VecOfVecXd us; // controls from last trajectory
  VecOfVecXd ls;
  VecOfMatXd Ls;
  double cost_s;

  // n = dims(state), m = dims(control)
  MatrixXd du; 	//m*T
  VecOfMatXd fx; //nxnx(T+1)
  VecOfMatXd fu; //nxmx(T+1)
  VecOfVecXd cx; //nx(T+1)
  VecOfVecXd cu; //mx(T+1)
  VecOfMatXd cxx; //nxnx(T+1)
  VecOfMatXd cxu; //nxmx(T+1)
  VecOfMatXd cuu; //mxmx(T+1)

  Vector2d dV; //2x1
	VecOfVecXd Vx; //nx(T+1)
	VecOfMatXd Vxx; //nxnx(T+1)
  VecOfVecXd k; //mxnxT
  VecOfMatXd K; //mxT

  double forward_pass(const VectorXd& x0, const VecOfVecXd& u);
  int backward_pass();
  VecOfVecXd add_bias_to_u(const VecOfVecXd &u, const VecOfVecXd &l, const double alpha);
  double get_gradient_norm(const VecOfVecXd& l, const VecOfVecXd& u);

  void compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);
  void get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u);

  FRIEND_TEST(ILQRSetup, dDynamicsTest);
  FRIEND_TEST(ILQRSetup, dCostTest);
  FRIEND_TEST(ILQRSetup, ddCostTest);
  FRIEND_TEST(ILQRSetup, ForwardPassTest);
};

#endif
