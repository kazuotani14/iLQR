#include "gtest/gtest.h"
#include <cmath>
#include "double_integrator.h"

using Eigen::VectorXd;
static const double eq_tol = 1e-6;

// class DoubleIntegratorSetup : public ::testing::Test
// {
// public:
//   Setup()
//   { }
//   DoubleIntegrator model;
// };

TEST(DoubleIntegratorTest, DxTest)
{
  DoubleIntegrator model;
  VectorXd x(4); x << 0.0, 0.0, 0.5, 0.1;
  VectorXd u(2); u << 1, -1;

  VectorXd dx = model.dynamics(x, u);
  VectorXd expected(4); expected << 0.5, 0.1, 1., -1.;
  EXPECT_TRUE(dx.isApprox(expected, eq_tol));

}

TEST(DoubleIntegratorTest, IntegrationTest)
{
  DoubleIntegrator model;
  VectorXd x(4); x << 0.0, 0.0, 0.5, 0.1;
  VectorXd u(2); u << 1, -1;
  double dt = 0.05;

  VectorXd dx = model.dynamics(x, u);
  VectorXd x1 = model.integrate_dynamics(x, u, dt);

  VectorXd expected = x + dt*dx;

  EXPECT_TRUE(x1.isApprox(expected, eq_tol));
}

TEST(DoubleIntegratorTest, CostTest)
{
  DoubleIntegrator model;
  VectorXd x(4); x << 0.1, 0.1, 0.5, 0.1;
  VectorXd u(2); u << 0.1, -1;

  double cost = model.cost(x, u);
  EXPECT_TRUE(std::abs(cost-2.682) < 0.001);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}