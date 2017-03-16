#include "loco_car.h"

int main(){
		std::cout << "started main\n";

		std::cout << "instantiating..\n";
	  LocoCar car;

		//   Eigen::VectorXd x(6);
		//   x << 1,2,3,1,2,3;
		//   Eigen::Vector2d u(1,1);
		//   Eigen::VectorXd dx(6);
		//   dx = loco.dynamics(x,u);
		//   print_vec(dx);
		//   std::cout << "Should be: \n -1.2722, -1.8389, 3.0000, 8.0671, -5.7746, 2.9479\n";

		Eigen::VectorXd x(8);
	  x << 1, 2, 3, 1, 2, 3, 0.1, 0.5;
	  Eigen::VectorXd u(2);
	  u << 1, 0.5;
	  Eigen::VectorXd x_d(6);
	  x_d << 3, 3, 0, 0, 0, 0;
	  Eigen::VectorXd obs(2);
	  obs << 3, 3;

		car.obs = obs;
		car.x_d = x_d;

		Eigen::VectorXd dx;
	  double c = car.get_dynamics_and_cost(x, u, dx);
	  std::cout << dx << '\n' << "cost: " << c;
}
