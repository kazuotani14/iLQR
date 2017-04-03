#include "ilqr.h"
#include "double_integrator.h"

int main()
{
	double dt = 0.05;
	int T = 50;
	VectorXd goal(4);
	goal << 1.0, 1.0, 0.0, 0.0;

	DoubleIntegrator* simple_model = new DoubleIntegrator(goal);

	iLQR ilqr_simple(simple_model, dt, T);

	VectorXd x0(4), xd(4);
	x0 << 0., 0., 0., 0.;

	VecOfVecXd u0;
	Vector2d u_init(0.1, 0.1);
	for (int i=0; i<T; i++)  u0.push_back(u_init);

	ilqr_simple.init_traj(x0,u0);
	ilqr_simple.generate_trajectory();
	std::cout << "Keep going!" << std::endl;

	return 0;
}
