#include "iLQR2.h"

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
												VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost)
{
	// Initialize dummy vectors;
	VecOfVecXd x;
	VecOfMatXd L;
	VecOfVecXd du;
	double alpha = 9999; //just so bugs are easier to identify

	forward_pass(x0, u, xnew, unew, new_cost, x, L, du, alpha);
}

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
									const VecOfVecXd &x, const VecOfMatXd &L,
									const VecOfVecXd &du, const double &alpha)
{

	double total_cost = 0;

	VecXd x_curr = x0;
	VecXd x1;
	VecXd u_curr;
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
			VecXd dx;
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
}
