#include "ilqr.h"

// initial roll-out
double iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u)
{
	// Initialize dummy vectors;
	VecOfVecXd x;
	VecOfMatXd L;

	return forward_pass(x0, u, x, L);
}

double iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u,
						const VecOfVecXd &x, const VecOfMatXd &L)
{
	double total_cost = 0;

	VectorXd x_curr = x0;
	VectorXd x1;
	VectorXd u_curr;
	xs[0] = x0;

	for(int t=0; t<T; t++)
	{
		u_curr = u[t];

		if (x.size()>0 && L.size()>0)
		{
			VectorXd dx = xs[t] - x[t];
			u_curr += L[t]*dx; //apply LQR control gains
		}

		us[t] = clamp_to_limits(u_curr, model->u_min, model->u_max);
		x1 = model->integrate_dynamics(x_curr, u_curr, dt);
		total_cost += model->cost(x_curr,u_curr);

		xs[t+1] = x1;
		x_curr = x1;
	}

	// calculate final cost
	total_cost += model->final_cost(xs[xs.size()-1]);

	return total_cost;
}
