#include "iLQR.h"

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
												VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost)
{
	// Initialize dummy vectors;
	VecOfVecXd x;

	forward_pass(x0, u, xnew, unew, new_cost, x;
}

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
									const VecOfVecXd &x)
{
	double total_cost = 0;

	VecXd x_curr = x0;
	VecXd x1;
	VecXd u_curr;
	xnew[0] = x0;

	for (int t=0; t<T; t++) 	//at each timestep
	{
		u_curr = u0[t];

		//clamp to min and max values in control_limits
		unew[t] = clamp_to_limits(u_curr);

		double cost = get_nextstate_and_cost(x_curr, u_curr, x1); //step forward in time
		total_cost += cost;

		xnew[t+1] = x1;
		x_curr = x1;
		// std::cout << x << "\n-----\n";
	}

	new_cost = total_cost;
}
