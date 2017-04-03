#include "ilqr.h"

double iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u)
{
	double total_cost = 0;

	VectorXd x_curr = x0;
	VectorXd x1;
	VectorXd u_curr;

	VecOfVecXd x_new(T+1);
	x_new[0] = x0;

	// for(const auto& ui : u) print_eigen("u_fp", ui);

	for(int t=0; t<T; t++)
	{
		u_curr = u[t];

		if (K.size()>0)
		{
			VectorXd dx = x_new[t] - xs[t];
			u_curr += K[t]*dx; //apply LQR control gains
		}

		us[t] = clamp_to_limits(u_curr, model->u_min, model->u_max);
		x1 = model->integrate_dynamics(x_curr, u_curr, dt);
		total_cost += model->cost(x_curr,u_curr);

		x_new[t+1] = x1;
		x_curr = x1;
	}

	// calculate final cost
	xs = x_new;
	total_cost += model->final_cost(xs[T]);
	// cout << model->final_cost(xs[T]) << endl;

	return total_cost;
}
