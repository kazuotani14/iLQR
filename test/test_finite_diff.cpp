#include "gtest/gtest.h"

#include "finite_diff.h"
#include "eigen/Eigen/Core"
#include <functional>
#include <iostream>
#include <cmath>

static const double eq_tol = 1e-6;

using Eigen::Vector2d;

TEST(FiniteDiffTest, Scalar2Scaler)
{
  std::function<double(double)> neg_quad = [](double x){return -pow(x,2);};
  double dx = finite_diff_scalar2scalar(neg_quad, -1);
  EXPECT_TRUE( (dx-2)<eq_tol );
}

TEST(FiniteDiffTest, Vec2Scalar)
{
  std::function<double(VectorXd)> quad_vec = [](VectorXd v){return v.array().square().sum();};
  VectorXd dx = finite_diff_vec2scalar(quad_vec, Vector2d(1.,1.));
  EXPECT_TRUE(dx.isApprox(Vector2d(2.,2.), eq_tol));
}

TEST(FiniteDiffTest, Vec2Vec)
{
  std::function<VectorXd(VectorXd)> add_ones = [](VectorXd v){return v + VectorXd::Ones(v.size()); };
  MatrixXd dx = finite_diff_vec2vec(add_ones, Vector2d(1.,0.), 2);

  EXPECT_TRUE(dx.isApprox(MatrixXd::Identity(2,2), eq_tol));
}

TEST(FiniteDiffTest, LargeProblem)
{
  std::function<VectorXd(VectorXd)> f = [](VectorXd v){return (v.array()*3).matrix();};
  VectorXd vec(1000);
  vec.setOnes();
  MatrixXd dx = finite_diff_vec2vec(f, vec, 1000);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
