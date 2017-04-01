#ifndef _ILQR_H_
#define _ILQR_H_

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
  int T;          // number of state transitions (#timesteps-1)
  double dt;

  // Tracking progress
  VecOfVecXd xs; // states from last trajectory
  VecOfVecXd us; // controls from last trajectory
  VecOfVecXd ls;
  VecOfMatXd Ls;

  // Stuff relevant at each cycle
  VecXd x_current;
  VecXd u_current;
  double cost_s;

  // Helper functions will modify passed-in variables instead of returning values,
  //  because most of them need to return multiple values.
  void forward_pass(const VecXd& x0, const VecOfVecXd& u,
  									VecOfVecXd& xnew, VecOfVecXd& unew, double& new_cost);

  void forward_pass(const VecXd& x0, const VecOfVecXd& u,
  									VecOfVecXd& xnew, VecOfVecXd& unew, double& new_cost,
  									const VecOfVecXd& x, const VecOfMatXd& L);

  int backward_pass(const VecOfVecXd &cx, const VecOfVecXd &cu, const VecOfMatXd &cxx, const VecOfMatXd &cxu,
  									const VecOfMatXd &cuu, const VecOfMatXd &fx, const VecOfMatXd &fu, const VecOfVecXd &u,
  									VecOfVecXd &Vx, VecOfMatXd &Vxx, VecOfVecXd &k, VecOfMatXd &K, Vec2d &dV);
  VecOfVecXd adjust_u(VecOfVecXd &u, VecOfVecXd &l, double alpha);

  // Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
  void compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd &fx,
  												 VecOfMatXd &fu, VecOfVecXd &cx, VecOfVecXd &cu,
  												 VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu);
  void get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
                                          VecOfMatXd &fx, VecOfMatXd &fu);
  void get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
  															   VecOfVecXd &cx, VecOfVecXd &cu);
  void get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
                                VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu);

public:
  iLQR(Model* p_dyn, double timeDelta, int timesteps): dt(timeDelta), T(timesteps-1)
  {
    model.reset(p_dyn);
  }

  VecXd x_d; // target state
  std::shared_ptr<Model> model;

  double init_traj(VecXd &x_0, VecOfVecXd &u0);
  void generate_trajectory();

  // Tester functions
  void output_to_csv();

};

#endif
