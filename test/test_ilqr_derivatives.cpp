#include "gtest/gtest.h"
#include "ilqr.h"
#include "double_integrator.h"
#include "eigen_helpers.h"

static const double eq_tol = 1e-2;

class ILQRSetup : public ::testing::Test
{
public:
  virtual void SetUp()
  {
    double dt = 0.05;
    int T = 10;
    VectorXd goal(4);
    goal << 1.0, 1.0, 0.0, 0.0;
    DoubleIntegrator* simple_model = new DoubleIntegrator(goal);

    ilqr_simple.reset(new iLQR(simple_model, dt));

    x0.resize(4);
    xd.resize(4);
    x0 << 0., 0., 0., 0.;

    Vector2d u_init(0.1, 0.1);
    for (int i=0; i<T; i++)  u0.push_back(u_init);

    ilqr_simple->init_traj(x0,u0);
  }

  VectorXd x0, xd;
  VecOfVecXd u0;

  std::shared_ptr<iLQR> ilqr_simple;
};


TEST_F(ILQRSetup, dDynamicsTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);
  Eigen::Matrix4d fx_expected;
  fx_expected << 1., 0., 0.05, 0.,
           0., 1., 0., 0.05,
           0., 0., 1., 0.,
           0., 0., 0., 1.;
  Eigen::MatrixXd fu_expected(4,2);
  fu_expected << 0., 0.,
            0., 0.,
            0.05, 0.,
            0., 0.05;
  // print_eigen("fx[0]", ilqr_simple->fx[0]);
  // print_eigen("fx expected", fx_expected);
  // print_eigen("fu[0]", ilqr_simple->fu[0]);
  // print_eigen("fu expected", fu_expected);
  EXPECT_TRUE(ilqr_simple->fx[0].isApprox(fx_expected, eq_tol));
  EXPECT_TRUE(ilqr_simple->fu[0].isApprox(fu_expected, eq_tol));
}

TEST_F(ILQRSetup, dCostTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);
  Eigen::Vector4d cx_expected;
  cx_expected << -2.0, -2.0, 0., 0.;
  Eigen::Vector2d cu_expected;
  cu_expected << 0.2, 0.2;
  // print_eigen("cx[3]", ilqr_simple->cx[3]);
  // print_eigen("cu[0]", ilqr_simple->cu[0]);
  EXPECT_TRUE(ilqr_simple->cx[0].isApprox(cx_expected, eq_tol));
  EXPECT_TRUE(ilqr_simple->cu[0].isApprox(cu_expected, eq_tol));
}

TEST_F(ILQRSetup, ddCostTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);

  Eigen::Matrix4d cxx_expected;
  cxx_expected << 2., 0.,0., 0.,
         0., 2., 0., 0.,
         0., 0., 0.4, 0.,
         0., 0., 0, 0.4;
  Eigen::MatrixXd cxu_expected(4,2);
  cxu_expected.setZero();
  Eigen::Matrix2d cuu_expected;
  cuu_expected << 2, 0,
         0, 2;

  // print_eigen("cxx[end]", ilqr_simple->cxx[ilqr_simple->cxx.size()-1]);
  // print_eigen("cxu[end]", ilqr_simple->cxu[ilqr_simple->cxu.size()-1]);
  // print_eigen("cuu[end]", ilqr_simple->cuu[ilqr_simple->cuu.size()-1]);
  // std::cout << ilqr_simple->cxx.size() << std::endl;
  EXPECT_TRUE(ilqr_simple->cxx[0].isApprox(cxx_expected, eq_tol));
  EXPECT_TRUE(ilqr_simple->cxu[0].isApprox(cxu_expected, eq_tol));
  EXPECT_TRUE(ilqr_simple->cuu[0].isApprox(cuu_expected, eq_tol));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
