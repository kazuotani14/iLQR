#ifndef _ILQR_H_
#define _ILQR_H_

#include "common.h"
#include "model.h"
#include "boxqp.h"

#include <memory>
#include <numeric>
#include <stdexcept>

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
  iLQR(Model* p_dyn, double timeDelta): dt(timeDelta)
  {
    model.reset(p_dyn);
  }
  iLQR() = default;

  std::shared_ptr<Model> model;

  void generate_trajectory();
  void generate_trajectory(const VectorXd &x_0); //warm-start
  void generate_trajectory(const VectorXd &x_0, const VecOfVecXd &u0); //fresh start

  void output_to_csv(const std::string filename);
  double init_traj(const VectorXd &x_0, const VecOfVecXd &u_0);

private:
  double dt;
  int T;  // number of state transitions

  VectorXd x0;

  VecOfVecXd xs; // s: "step". current working sequence
  VecOfVecXd us;
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
  void get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& f_x, VecOfMatXd& f_u);
  void get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfVecXd& c_x, VecOfVecXd& c_u);
  void get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& c_xx, VecOfMatXd& c_xu, VecOfMatXd& c_uu);

  //multi-threaded version - about 2 to 3 times faster
  void get_cost_2nd_derivatives_mt(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& c_xx, VecOfMatXd& c_xu, VecOfMatXd& c_uu, int n_threads_per);

  void calculate_cxx(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& c_xx, int start_T, int end_T);
  void calculate_cxu(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& c_xu, int start_T, int end_T);
  void calculate_cuu(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd& c_uu, int start_T, int end_T);

  FRIEND_TEST(ILQRSetup, dDynamicsTest);
  FRIEND_TEST(ILQRSetup, dCostTest);
  FRIEND_TEST(ILQRSetup, ddCostTest);
  FRIEND_TEST(ILQRSetup, ForwardPassTest);

};

#endif
