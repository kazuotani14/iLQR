#include "iLQR2.h"

double iLQR::init_traj(VecXd &x_0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.reserve(T+1);
	us.reserve(T);

	// call forward_pass, get xs, us, cost
	double initial_cost;
	forward_pass(x_0, u0, xs, us, initial_cost);
	std::cout << "Initial cost: " << initial_cost << "\n";

	return initial_cost;
}


void iLQR::generate_trajectory(const VecXd &x_0, const VecXd &x_d, const int trajectoryLength)
{
	// TODO Check inputs

	// TODO Make sure trajectory (xs, us) is initialized, copy to x and u
	VecOfVecXd x;	//nxT
	VecOfVecXd u; //2xT // TODO change these to vectors

	// Initialize all vectors, matrices we'll be using
	VecOfMatXd du(T); //2*T double
	VecOfMatXd fx(T+1); //nxnx(T+1)
	VecOfMatXd fu(T+1); //nx2x(T+1)
	VecOfVecXd cx(T+1); //nx(T+1)
	VecOfVecXd cu(T+1); //2x(T+1)
	VecOfMatXd cxx(T+1); //nxnx(T+1)
	VecOfMatXd cxu(T+1); //nx2x(T+1)
	VecOfMatXd cuu(T+1); //2x2x(T+1), nxnx(T+1)
	//fxx: none, fxu: None, fuu: none (because we're only doing first-order)

	VecOfVecXd Vx; //nx(T+1)
	VecOfMatXd Vxx; //nxnx(T+1)
	VecOfMatXd L(T); //2xnxT
	VecOfVecXd l(T); //2xT
	Vec2d dV; //2x1

	// constants, timers, counters
	bool flgChange = true;
	bool stop = false;
	double dcost = 0;
	double z = 0;
	double expected = 0;
	std::cout << "\n=========== begin iLQG ===========\n";

	for (int iter=1; iter<maxIter; iter++)
	{
		if (stop)
			break;

		//--------------------------------------------------------------------------
		// STEP 1: Differentiate dynamics and cost along new trajectory
		if (flgChange){
			compute_derivatives(xs,us, fx,fu,cx,cu,cxx,cxu,cuu);
			flgChange = 0;
		}

		//--------------------------------------------------------------------------
		// STEP 2: Backward pass, compute optimal control law and cost-to-go
		bool backPassDone = false;
		while (!backPassDone)
		{
	 		// update Vx, Vxx, l, L, dV with back_pass
			// TODO
			int diverge = backward_pass(cx,cu,cxx,cxu,cuu,fx,fu,u, Vx, Vxx, l, L, dV);

		//
		// 	if (diverge!=0)
		// 	{
		// 		std::cout << "Cholesky failed at timestep " << diverge << ".\n";
		// 		dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
		// 		lambda    = std::max(lambda * dlambda, lambdaMin);
		// 		if lambda > lambdaMax
		// 				break;
		// 		continue;
		// 	}
			backPassDone = true;
		}
		// TODO check for termination due to small gradient

		//--------------------------------------------------------------------------
		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost
		// bool fwdPassDone = 0;
		// if (backPassDone)
		// {
		// 	TODO this whole section
		// 	TODO initialize xnew, unew, costnew
		// 	[xnew,unew,costnew] = forward_pass(x0 ,u, L, x(:,1:N), l, Op.Alpha, DYNCST,Op.lims,Op.diffFn);
		// 	Dcost               = sum(cost(:)) - sum(costnew,2);
		// 	[dcost, w]          = max(Dcost); % should be positive
		// 	alpha               = Alpha(w);
		// 	expected            = -alpha*(dV(1) + alpha*dV(2));
		// 	if (expected > 0){
		// 		z = dcost/expected;
		// 	}
		// 	else{
		// 		z = sign(dcost);
		// 		std::cout << "non-positive expected reduction: should not occur\n";
		// 	}
		// 	if (z > Op.zMin){
		// 		fwdPassDone = 1;
		// 		costnew     = costnew(:,:,w);
		// 		xnew        = xnew(:,:,w);
		// 		unew        = unew(:,:,w);
		// 	}
		// }
		// if (!fwdPassDone){
		// 	alpha = NULL; // signals failure of forward pass
		// }

	// 	//--------------------------------------------------------------------------
	// 	// STEP 4: accept step (or not), print status
	// 	if (iter==1)
	// 		std::cout << "iteration\tcost\treduction\texpected\tgradient\tlog10(lambda)\n";
	//
	// 	if (fwdPassDone)
	// 	{
	// 		// TODO print stuff
	//
	// 		// decrease lambda
	// 		dlambda   = min(dlambda / lambdaFactor, 1/lambdaFactor);
	// 		lambda    = lambda * dlambda * (lambda > lambdaMin);
	//
	// 		// accept changes
	// 		u              = unew;
	// 		x              = xnew;
	// 		cost           = costnew;
	// 		flgChange      = true;
	//
	// 		// terminate?
	// 		if (dcost < tolFun){
	// 			std::cout << "\nSUCCESS: cost change < tolFun\n";
	// 			break;
	// 		}
	// 	}
	// 	else // no cost improvement
	// 	{
	// 		// increase lambda
	// 		dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
	// 		lambda   = std::max(lambda * dlambda, lambdaMin);
	//
	// 		// TODO print stuff
	//
	// 		// terminate?
	// 		if (lambda > lambdaMax){
	// 			std::cout << "\nEXIT: lambda > lambdaMax\n";
	// 			break;
	// 		}
	// 	}
	} // optimization for loop
	//
	// if (iter==maxIter)
	// 	std::cout << "\nEXIT: Maximum iterations reached.\n";
	//
	// // TODO print final stuff

}//generate_trajectory
