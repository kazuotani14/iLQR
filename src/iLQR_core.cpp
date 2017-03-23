#include "iLQR.h"

#define VERBOSE 1

double iLQR::init_traj(VecXd &x_0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.reserve(T+1);
	us.reserve(T);

	// call forward_pass, get xs, us, cost
	double cost_i;
	forward_pass(x_0, u0, xs, us, cost_i);
	std::cout << "Initial cost: " << cost_i << "\n";

	initial_cost = cost_i;

	return initial_cost;
}

// double iLQR::get_gradient_norm(VecOfVecXd l, VecOfVecXd u)
// {
// 	for (int i=0; i<l.size())
// }

void iLQR::generate_trajectory(const VecXd &x_0, const VecXd &x_d, const int trajectoryLength)
{
	// TODO Check inputs

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
		std::cout << "Finished step 1 : compute derivatives. \n";

		//--------------------------------------------------------------------------
		// STEP 2: Backward pass, compute optimal control law and cost-to-go
		bool backPassDone = false;
		while (!backPassDone)
		{
	 		// update Vx, Vxx, l, L, dV with back_pass
			// TODl
			int diverge = backward_pass(cx,cu,cxx,cxu,cuu,fx,fu,u, Vx,Vxx,l,L,dV);
			// std::cout << "l: \n" << l[0] << '\n';
			// std::cout << "L: \n" << L[0] << '\n';

			if (diverge!=0)
			{
				std::cout << "Cholesky failed at timestep " << diverge << ".\n";
				dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
				lambda    = std::max(lambda * dlambda, lambdaMin);
				if (lambda > lambdaMax)
						break;
				continue;
			}
			backPassDone = true;
		}
		// TODO check for termination due to small gradient
		// double gnorm = get_gradient_norm(l, u);

		std::cout << "Finished step 2 : backward pass. \n";

		//--------------------------------------------------------------------------
		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost
		// 	TODO this whole section
		bool fwdPassDone = 0;
		// 	TODO initialize xnew, unew, costnew
		VecOfVecXd xnew(trajectoryLength+1);
		VecOfVecXd unew(trajectoryLength);
		double new_cost;

		if (backPassDone) //  serial backtracking line-search
		{
			for (int i=0; i<Alpha.size(); i++){
				forward_pass(x_0, u, xnew,unew,new_cost, x,L);
				dcost    = initial_cost - new_cost;
			}

		}
		// expected = -alpha*(dV(1) + alpha*dV(2));
		// if expected > 0
		// 		z = dcost/expected;
		// else
		// 		z = sign(dcost);
		// 		warning('non-positive expected reduction: should not occur');
		// end
		// if (z > Op.zMin)
		// 		fwdPassDone = 1;
		// 		break;

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
