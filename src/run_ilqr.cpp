#include "common.h"
#include "ilqr.h"
#include "double_integrator.h"
#include "acrobot.h"

int main(int argc, char *argv[]) {

  if (argc != 2 || 
      (strcmp(argv[1], "acrobot")!=0 && strcmp(argv[1], "integrator")!=0))  {
    cout << "Provide command line argument 'acrobot' or 'integrator'" << endl;
    return 0;
  }

  iLQR* ilqr;
  VecOfVecXd u0;
  VectorXd x0;

  if(strcmp(argv[1], "integrator") == 0) {
    // Make double integrator with goal state
    VectorXd goal(4);
    goal << 1.0, 0.5, 0.0, 0.0;
    Model* particle = new DoubleIntegrator(goal);
    // TODO make this part of iLQR problem rather than "dynamics model"?

    // Define problem
    double dt = 0.02;
    ilqr = new iLQR(particle, dt);

    // Define initial state
    x0.resize(4);
    x0 << -1.0, 0., 0.0, -0.2;

    // Make initialization for control sequence
    int T = 99; // right now, (T+1) has to be divisible by 10 - see derivatives.cpp. TODO remove this constraint
    Vector2d u_init; u_init.setZero();
    for (int i=0; i<T; i++) u0.push_back(u_init);
  }

  if(strcmp(argv[1], "acrobot") == 0) {
    // Make acrobot
    Model* acrobot = new Acrobot();

    // Define problem
    double dt = 0.02; // needs to be small enough for euler integration to be stable
    ilqr = new iLQR(acrobot, dt);

    // Define initial state
    x0.resize(4); x0.setZero();

    // Make initialization for control sequence
    int T = 499; 
    Vector1d u_init; u_init.setZero();
    for (int i=0; i<T; i++) u0.push_back(u_init);
  }

  // Solve!
  cout << "Run iLQR!" << endl;
  auto start = std::chrono::system_clock::now();
  ilqr->generate_trajectory(x0, u0);
  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  cout << "iLQR took: " << elapsed/1000. << " seconds." << endl;

  return 0;
}
