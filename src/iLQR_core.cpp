#include "iLQR.h"

#define VERBOSE 1

double iLQR::init_traj(VecXd &x_0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.resize(T+1);
	us.resize(T);

	// call forward_pass, get xs, us, cost
	double cost_i;
	forward_pass(x_0, u0, xs, us, cost_i);
	std::cout << "Initial cost: " << cost_i << "\n";

	cost_s = cost_i;

	return cost_s;
}


VecOfVecXd iLQR::adjust_u(VecOfVecXd &u, VecOfVecXd &l, double alpha)
{
	VecOfVecXd new_u = u;
	for(int i=0; i<u.size(); i++){
		new_u[i] += l[i]*alpha;
	}
	return new_u;
}

// double iLQR::get_gradient_norm(VecOfVecXd l, VecOfVecXd u)
// {
// 	for (int i=0; i<l.size())
// }

void iLQR::generate_trajectory(const VecXd &x_0, const VecXd &x_d, const int trajectoryLength)
{
	// TODO Check inputs
	if (us.size()==0){
		std::cout << "Call init_traj first.\n";
		return;
	}

	// TODO Make sure trajectory (xs, us) is initialized, copy to x and u
	VecOfVecXd x;	//nxT
	VecOfVecXd u; //2xT // TODO change these to vectors

	// Initialize all vectors, matrices we'll be using
	MatXd du(2,T); 	//2*T double
	VecOfMatXd fx(T+1); //nxnx(T+1)
	VecOfMatXd fu(T+1); //nx2x(T+1)
	VecOfVecXd cx(T+1); //nx(T+1)
	VecOfVecXd cu(T+1); //2x(T+1)
	VecOfMatXd cxx(T+1); //nxnx(T+1)
	VecOfMatXd cxu(T+1); //nx2x(T+1)
	VecOfMatXd cuu(T+1); //2x2x(T+1)

	VecOfVecXd Vx(T+1); //nx(T+1)
	VecOfMatXd Vxx(T+1); //nxnx(T+1)
	VecOfMatXd L(T); //2xnxT
	VecOfVecXd l(T); //2xT
	Vec2d dV; //2x1

	for (int i=0; i<=T+1; i++){
		fx[i].resize(n,n);
		fu[i].resize(n,m);
		cx[i].resize(n);
		cu[i].resize(m);
		cxx[i].resize(n,n);
		cxu[i].resize(n,m);
		cuu[i].resize(m,m);
		Vx[i].resize(n);
		Vxx[i].resize(n,n);
		if(i<=T){
			L[i].resize(2,n);
			l[i].resize(2);
		}
	}

	// constants, timers, counters
	bool flgChange = true;
	bool stop = false;
	double dcost = 0;
	double z = 0;
	double expected = 0;
	int diverge = 0;
	std::cout << "\n=========== begin iLQG ===========\n";

	int iter;
	for (iter=0; iter<maxIter; iter++)
	{
		if (stop)
			break;

		// std::cout << "Iteration " << iter << ".\n";
		//--------------------------------------------------------------------------
		// STEP 1: Differentiate dynamics and cost along new trajectory
		if (flgChange){
			compute_derivatives(xs,us, fx,fu,cx,cu,cxx,cxu,cuu);
			flgChange = 0;
		}
		// std::cout << "Finished step 1 : compute derivatives. \n";

		//--------------------------------------------------------------------------
		// STEP 2: Backward pass, compute optimal control law and cost-to-go
		bool backPassDone = false;
		while (!backPassDone)
		{
	 		// update Vx, Vxx, l, L, dV with back_pass
			diverge = 0;
			diverge = backward_pass(cx,cu,cxx,cxu,cuu,fx,fu,u, Vx,Vxx,l,L,dV);
			// std::cout << "l: \n" << l[0] << '\n';
			// std::cout << "L: \n" << L[0] << '\n';

			if (diverge!=0)
			{
				// std::cout << "Cholesky failed at timestep " << diverge << ".\n";
				dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
				lambda    = std::max(lambda * dlambda, lambdaMin);
				if (lambda > lambdaMax)
						break;
				continue;
			}
			backPassDone = true;
		}
		// TODO check for termination due to small gradient
		double gnorm = 0;
		// double gnorm = get_gradient_norm(l, u);

		// std::cout << "Finished step 2 : backward pass. \n";

		//--------------------------------------------------------------------------
		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost
		bool fwdPassDone = 0;
		VecOfVecXd xnew(trajectoryLength+1);
		VecOfVecXd unew(trajectoryLength);
		double new_cost, alpha;

		if (backPassDone) //  serial backtracking line-search
		{
			for (int i=0; i<Alpha.size(); i++){
				alpha = Alpha(i);
				forward_pass(x_0, adjust_u(us,l,alpha), xnew,unew,new_cost, x,L);
				dcost    = cost_s - new_cost;
				expected = -alpha * (dV(0) + alpha*dV(1));

				// std::cout << "dV: " << dV(0) << ' ' << dV(1) << '\n';
				// std::cout << "dcost: " << dcost << '\n';
				// std::cout << "expected: " << expected << '\n';
				// getchar();

				if (expected>0){
					z = dcost/expected;
				}
				else{
					z = sgn(dcost);
					std::cout << "Warning: non-positive expected reduction: should not occur\n";
				}
				if(z > zMin){
					fwdPassDone = 1;
					break;
				}
			}

			if (!fwdPassDone){
				alpha = 0.0; // signals failure of forward pass
			}
		}

		// std::cout << "Finished step 3 : forward pass. \n";

	//--------------------------------------------------------------------------
	// STEP 4: accept step (or not), print status
	 	if (iter==0)
	 		std::cout << "iteration\tcost\t\treduction\texpected\tgradient\tlog10(lambda)\n";

	 	if (fwdPassDone)
	 	{
	 		// TODO print stuff
			printf("%-12d\t%-12.6g\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                iter, new_cost, dcost, expected, gnorm, log10(lambda));

	 		// decrease lambda
	 		dlambda   = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
	 		lambda    = lambda * dlambda * (lambda > lambdaMin);

	 		// accept changes
	 		us              = unew;
	 		xs              = xnew;
	 		cost_s          = new_cost;
	 		flgChange       = true;

	 		//terminate?
	 		if (dcost < tolFun){
	 			std::cout << "\nSUCCESS: cost change < tolFun\n";
	 			break;
	 		}
	 	}
	 	else // no cost improvement
	 	{
	 		// increase lambda
	 		dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
	 		lambda   = std::max(lambda * dlambda, lambdaMin);

	 		// TODO print stuff
			printf("%-12d\t%-12s\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                iter,"NO STEP", dcost, expected, gnorm, log10(lambda));
	 		// terminate?
	 		if (lambda > lambdaMax){
	 			std::cout << "\nEXIT: lambda > lambdaMax\n";
	 			break;
	 		}
		}
	} // optimization for loop
	//
	if (iter==maxIter)
	 		std::cout << "\nEXIT: Maximum iterations reached.\n";
	//
	// TODO print final stuff

} //generate_trajectory
