#include "iLQR2.h"

void iLQR::forward_pass(Eigen::VectorXd &x0, std::vector<Eigen::Vector2d> &u,
									std::vector<Eigen::VectorXd> &xnew, std::vector<Eigen::Vector2d> &unew, double &new_cost)
{
	// Initialize zero vectors;
	std::vector<Eigen::VectorXd> x;
	std::vector<Eigen::MatrixXd> L;
	std::vector<Eigen::Vector2d> du;
	double alpha = 9999; //just so bugs are easier to identify

	forward_pass(x0, u, xnew, unew, new_cost, x, L, du, alpha);
}

void iLQR::forward_pass(Eigen::VectorXd &x0, std::vector<Eigen::Vector2d> &u,
									std::vector<Eigen::VectorXd> &xnew, std::vector<Eigen::Vector2d> &unew, double &new_cost,
									std::vector<Eigen::VectorXd> &x, std::vector<Eigen::MatrixXd> &L,
									std::vector<Eigen::Vector2d> &du, double &alpha)
{
/*
INPUTS
 x0: nx1          u: 2xT            L: 2xnxT          x: nxT
 du: 2*T 	alpha: 1x1 or 11x1
OUTPUTS
 xnew, unew, costnew
*/

	double total_cost = 0;

	Eigen::VectorXd x_curr = x0;
	Eigen::VectorXd x1;
	Eigen::Vector2d u_curr;
	xnew[0] = x0;

	for (int t=0; t<T; t++) 	//at each timestep
	{
		u_curr = u0[t];

		if (du.size()>0)
		{
			u_curr += du[t]*alpha; // change u by du*alpha
		}

		if (L.size()>0 && x.size()>0)
		{
			Eigen::VectorXd dx;
			dx = xnew[t] - x[t];
			u_curr += L[t]*dx; //apply LQR control gains
		}

		//TODO clamp to min and max values in control_limits
		unew[t] = u_curr;

		double cost = get_dynamics_and_cost(x_curr, u_curr, x1); //step forward in time
		total_cost += cost;

		xnew[t+1] = x1;
		x_curr = x1;
		// std::cout << x << "\n-----\n";
	}

	new_cost = total_cost;
} //forward_pass
