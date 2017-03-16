#include "iLQR2.h"

void iLQR::init_traj(Eigen::VectorXd &x_0, std::vector<Eigen::Vector2d> &u0)
{
  // std::cout << "starting init_traj\n";
	//initialize xs and us, return or store initial cost
  double total_cost = 0;

  Eigen::VectorXd x = x_0;
  Eigen::VectorXd x1;
  Eigen::Vector2d u;
  xs.push_back(x_0);

  // std::cout << "entering loop\n";
	for (int t=0; t<T; t++) 	//at each timestep
  {
    u = u0[t];
    us.push_back(u);
    double cost = get_dynamics_and_cost(x, u, x1);
    total_cost += cost;
    xs.push_back(x1);
    x = x1;
    // std::cout << x << "\n-----\n";
  }
  // std::cout << "done with init_traj\n";
}
