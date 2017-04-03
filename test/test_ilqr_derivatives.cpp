#include "gtest/gtest.h"
#include "ilqr.h"
#include "double_integrator.h"
#include "eigen_helpers.h"

class ILQRSetup : public ::testing::Test
{
public:
  virtual void SetUp()
  {
    double dt = 0.05;
    int T = 50;
    VectorXd goal(4);
    goal << 1.0, 1.0, 0.0, 0.0;
    DoubleIntegrator* simple_model = new DoubleIntegrator(goal);

    ilqr_simple.reset(new iLQR(simple_model, dt, T));

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

//TODO add assertions

TEST_F(ILQRSetup, dDynamicsTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);
  print_eigen("fx[0]", ilqr_simple->fx[0]);
  print_eigen("fu[0]", ilqr_simple->fu[0]);
}

TEST_F(ILQRSetup, dCostTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);
  print_eigen("cx[0]", ilqr_simple->cx[0]);
  print_eigen("cu[0]", ilqr_simple->cu[0]);
}

TEST_F(ILQRSetup, ddCostTest)
{
  ilqr_simple->compute_derivatives(ilqr_simple->xs, ilqr_simple->us);
  print_eigen("cxx[0]", ilqr_simple->cxx[0]);
  print_eigen("cxu[0]", ilqr_simple->cxu[0]);
  print_eigen("cuu[0]", ilqr_simple->cuu[0]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
