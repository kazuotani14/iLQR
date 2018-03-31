#include "common.h"
#include "ilqr.h"
#include "double_integrator.h"
#include "acrobot.h"

// TODO do this with command line arguments
// #define DOUBLEINTEGRATOR
#define ACROBOT

int main() {

#ifdef DOUBLEINTEGRATOR
  // Make double integrator with goal state
  VectorXd goal(4);
  goal << 1.0, 0.5, 0.0, 0.0;
  Model* particle = new DoubleIntegrator(goal);
  // TODO make this part of iLQR problem rather than "dynamics model"?

  // Define problem
  double dt = 0.05;
  iLQR ilqr_simple(particle, dt);

  // Define initial state
  VectorXd x0(4);
  x0 << -1.0, 0., 0.0, -0.2;

  // Make initialization for control sequence
  int T = 99; // right now, (T+1) has to be divisible by 10 - see derivatives.cpp. TODO remove this constraint
  VecOfVecXd u0;
  Vector2d u_init; u_init.setZero();
  for (int i=0; i<T; i++) u0.push_back(u_init);
#endif

#ifdef ACROBOT
  // Make acrobot
  Model* acrobot = new Acrobot();

  // Define problem
  double dt = 0.02;
  iLQR ilqr_simple(acrobot, dt);

  // Define initial state
  VectorXd x0(4); x0.setZero();

  // Make initialization for control sequence
  int T = 499; 
  VecOfVecXd u0;
  Vector1d u_init; u_init.setZero();
  for (int i=0; i<T; i++) u0.push_back(u_init);
#endif

  // Solve!
  cout << "Run iLQR!" << endl;
  auto start = std::chrono::system_clock::now();
  ilqr_simple.generate_trajectory(x0, u0);
  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  cout << "iLQR took: " << elapsed/1000. << " seconds." << endl;

  return 0;
}
