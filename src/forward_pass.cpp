#include "ilqr.h"

// initial roll-out
void iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u,
						VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost)
{
	// Initialize dummy vectors;
	VecOfVecXd x;
	VecOfMatXd L;

	forward_pass(x0, u, xnew, unew, new_cost, x, L);
}

void iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u,
						VecOfVecXd &xnew, VecOfVecXd &unew, double &new_cost,
						const VecOfVecXd &x, const VecOfMatXd &L)
{
	double total_cost = 0;

	VectorXd x_curr = x0;
	VectorXd x1;
	VectorXd u_curr;
	xnew[0] = x0;

	for(int t=0; t<T; t++)
	{
		u_curr = u[t];

		if (x.size()>0 && L.size()>0)
		{
			VectorXd dx = xnew[t] - x[t];
			u_curr += L[t]*dx; //apply LQR control gains
		}

		//clamp to min and max values in control_limits
		unew[t] = clamp_to_limits(u_curr, model->u_min, model->u_max);

		x1 = model->integrate_dynamics(x_curr, u_curr, dt); //step forward in time
		total_cost += model->cost(x_curr,u_curr);

		xnew[t+1] = x1;
		x_curr = x1;
	}

	// calculate final cost
	total_cost += model->final_cost(xnew[xnew.size()-1]);

	new_cost = total_cost;
}
