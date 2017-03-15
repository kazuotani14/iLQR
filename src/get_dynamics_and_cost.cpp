#include "iLQR2.h"

void iLQR::finite_difference()
{
	// need function pointer
	// or just make separate functions for finite differencing different things, like TGlad did
}


void iLQR::car_dynamics_and_cost(Eigen::Matrix2D x, Eigen::Matrix2D u,

	Eigen::Matrix2D f, double c, Eigen::MatrixXD fx, Eigen::MatrixXD fu,
	Eigen::Matrix2D cx,Eigen::Matrix2D cu, Eigen::MatrixXD cxx,
	Eigen::MatrixXD cxu, Eigen::MatrixXD cuu, int mode)
{
	// This has a weird condition based on number of outputs in MATLAB code
	// 	What's the best way to deal with this? - Split into two functions
	if (mode==1)
	{
		f = dynamics(x,u);
		c = car_cost(x,u);
	}
	else if (mode==2)
	{
		//TODO pass in past control input - how? just make another member variable?


	}
	else
	{
		std::cout << "Select mode 1 or 2!\n";
	}
}
