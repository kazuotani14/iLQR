#include "loco_car.h"

int main(){
	  LocoCar car;

		VecXd x0(6);
		x0 << 0, 0, 0, 2, 0, 0;
		VecXd x_d(6);
		x_d << 3, 3, 0, 0, 0, 0;
		VecXd obs(2);
		obs << 3, 3;
		car.obs = obs;
		car.x_d = x_d;

		int T = 50;

		VecOfVecXd u0;
		Vec2d u_init(1, 0.3);
		for (int i=0; i<T; i++){
			u0.push_back(u_init);
		}

		car.u0 = u0;
		car.init_traj(x0,u0);
		// car.demoQP();
		car.generate_trajectory(x0, x_d, T);
}

//   VecXd x(6);
//   x << 1,2,3,1,2,3;
//   Vec2d u(1,1);
//   VecXd dx(6);
//   dx = loco.dynamics(x,u);
//   print_vec(dx);
//   std::cout << "Should be: \n -1.2722, -1.8389, 3.0000, 8.0671, -5.7746, 2.9479\n";

// VecXd x(8);
// x << 1, 2, 3, 1, 2, 3, 0.1, 0.5;
// VecXd u(2);
// u << 1, 0.5;


// VecXd x1;
// double c = car.get_nextstate_and_cost(x, u, x1);
// std::cout << x1 << '\n' << "cost: " << c << '\n';
