#include "iLQR.h"

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
												VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost)
{
	// Initialize dummy vectors;
	VecOfVecXd x;
	VecOfMatXd L;


	forward_pass(x0, u, xnew, unew, new_cost, x, L);
}

void iLQR::forward_pass(const VecXd &x0, const VecOfVecXd &u,
									VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
									const VecOfVecXd &x, const VecOfMatXd &L)
{
	double total_cost = 0;

	VecXd x_curr = x0;
	VecXd x1;
	VecXd u_curr;
	xnew[0] = x0;

	int t;
	for (t=0; t<T; t++) 	//at each timestep
	{
		u_curr = u[t];

		if (x.size()>0 && L.size()>0)
		{
		VecXd dx;
			dx = xnew[t] - x[t];
			u_curr += L[t]*dx; //apply LQR control gains
		}

		//clamp to min and max values in control_limits
		unew[t] = clamp_to_limits(u_curr);

		x1 = integrate_dynamics(x_curr,u_curr); //step forward in time
		total_cost += cost(x_curr,u_curr);

		xnew[t+1] = x1;
		x_curr = x1;
	}

	// calculate final cost
	total_cost += final_cost(xnew[t]);

	new_cost = total_cost;
}
