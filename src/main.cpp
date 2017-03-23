#include "loco_car.h"

int main(){
	  LocoCar car;

		VecXd x0(6);
		x0 << 0, 0, 0, 2, 0, 0;
		 //x0 << 0.1, 0, 0, 1.94, 0.05232, 0.559;
		VecXd x_d(6);
		x_d << 3, 0, 0, 0, 0, 0;
		VecXd obs(2);
		obs << 1, 0;
		car.obs = obs;
		car.x_d = x_d;

		int T = 50;

		VecOfVecXd u0;
		Vec2d u_init(1, 0.3);
		for (int i=0; i<T; i++){
			u0.push_back(u_init);
			// u0.push_back(Vec2d::Random());
		}
		car.u0 = u0;

		std::clock_t start;
		start = std::clock();

		car.init_traj(x0,u0);
		car.generate_trajectory(x0, x_d, T);

		double time_elapsed = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
		std::cout << "Took " << time_elapsed << " seconds.\n";
}


// std::clock_t start;
// start = std::clock();
// double time_elapsed = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
// std::cout << "Took " << time_elapsed << " seconds.\n";

// Dynamics check
// VecXd x(6);
// x << 1,2,3,1,2,3;
// Vec2d u(1,1);
// VecXd dx(6);
// dx = car.dynamics(x,u);
// print_vec(dx);
// std::cout << "Should be: \n    -1.2722 -1.8389 3.0000 7.2273 -4.6562 3.2690\n";

// VecXd x(8);
// x << 1, 2, 3, 1, 2, 3, 0.1, 0.5;
// VecXd u(2);
// u << 1, 0.5;
