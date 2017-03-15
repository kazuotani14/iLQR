/*
BIBTeX:
@INPROCEEDINGS{
author={Tassa, Y. and Mansard, N. and Todorov, E.},
booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
title={Control-Limited Differential Dynamic Programming},
year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}
*/

#include "iLQR2.h"

void iLQR::init_traj(Eigen::VectorXd &x_0, Eigen::Matrix2d &u0)
{
	std::cout << "Performing initial rollout.\n";
	// TODO initialize x_current and u_current

	// TODO roll out using dynamics
}


void iLQR::generate_trajectory(const Eigen::VectorXd &x_0, const Eigen::VectorXd &x_d,
																int trajectoryLength, int nControlInputs)
{
	/*
	 * x_0: start state
	 * x_d: goal state
	*/

	// TODO Check inputs

	// TODO Make sure trajectory (xs, us) is initialized, copy to x and u
	Eigen::Matrix2d x;	//nxT
	Eigen::Matrix2d u; //2xT

	// Initialize all vectors, matrices we'll be using
	Eigen::MatrixXd L; //2xnxT
	Eigen::Matrix2d l; //2xT
	Eigen::Matrix2d du; //2*T double
	Eigen::Vector2d alpha; //1x1 or 11x1
	Eigen::Matrix2d cx; //2x(T+1)
	Eigen::Matrix2d cu: //2x(T+1)
	Eigen::MatrixXd cuu; //nxnx(T+1)
	Eigen::MatrixXd cxu; //?
	Eigen::MatrixXd cxx; //nxnx(T+1)
	Eigen::MatrixXd cuu: //2x2x(T+1)
	Eigen::MatrixXd fx; //nxnx(T+1)
	Eigen::MatrixXd fu; //nx2x(T+1)
	//fxx: none, fxu: None, fuu: none (because we're only doing first-order)
	Eigen::Matrix2d Vx; //nx(T+1)
	Eigen::Matrix2d Vxx; //nxnx(T+1)
	Eigen::Vector2d dV; //2x1

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
		// STEP 1: Forward pass, differentiate dynamics and cost along new trajectory
		if flgChange{
			// TODO update fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu with DYNCST
			car_dynamics_and_cost(x,u,  f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu, 1);
			flgChange = 0;
		}

		//--------------------------------------------------------------------------
		// STEP 2: Backward pass, compute optimal control law and cost-to-go
		bool backPassDone = false;
		while (!backPassDone)
		{
			// update Vx, Vxx, l, L, dV with back_pass
			int diverge = backward_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,u); // TODO

			if (diverge!=0)
			{
				std::cout << "Cholesky failed at timestep " << diverge << ".\n";
				dlambda   = max(dlambda * lambdaFactor, lambdaFactor);
				lambda    = max(lambda * dlambda, lambdaMin);
				if lambda > lambdaMax
						break;
				continue;
			}
			backPassDone = true;
		}
		// TODO check for termination due to small gradient

		//--------------------------------------------------------------------------
		// STEP 3: line-search to find new control sequence, trajectory, cost
		bool fwdPassDone = 0;
		if (backPassDone)
		{
			//TODO this whole section
			// TODO initialize xnew, unew, costnew
			[xnew,unew,costnew] = forward_pass(x0 ,u, L, x(:,1:N), l, Op.Alpha, DYNCST,Op.lims,Op.diffFn);
			Dcost               = sum(cost(:)) - sum(costnew,2);
			[dcost, w]          = max(Dcost); % should be positive
			alpha               = Alpha(w);
			expected            = -alpha*(dV(1) + alpha*dV(2));
			if (expected > 0){
				z = dcost/expected;
			}
			else{
				z = sign(dcost);
				std::cout << "non-positive expected reduction: should not occur\n";
			}
			if (z > Op.zMin){
				fwdPassDone = 1;
				costnew     = costnew(:,:,w);
				xnew        = xnew(:,:,w);
				unew        = unew(:,:,w);
			}
		}
		if (!fwdPassDone){
			alpha = NULL; // signals failure of forward pass
		}

		//--------------------------------------------------------------------------
		// STEP 4: accept step (or not), print status
		if (iter==1)
			std::cout << "iteration\tcost\treduction\texpected\tgradient\tlog10(lambda)\n";

		if (fwdPassDone)
		{
			// TODO print stuff

			// decrease lambda
			dlambda   = min(dlambda / lambdaFactor, 1/lambdaFactor);
			lambda    = lambda * dlambda * (lambda > lambdaMin);

			// accept changes
			u              = unew;
			x              = xnew;
			cost           = costnew;
			flgChange      = true;

			// terminate?
			if (dcost < tolFun){
				std::cout << "\nSUCCESS: cost change < tolFun\n";
				break;
			}
		}
		else // no cost improvement
		{
			// increase lambda
			dlambda  = max(dlambda * lambdaFactor, lambdaFactor);
			lambda   = max(lambda * dlambda, lambdaMin);

			// TODO print stuff

			// terminate?
			if (lambda > lambdaMax){
				std::cout << "\nEXIT: lambda > lambdaMax\n";
				break;
			}
		}
	} // optimization for loop

	if (iter==maxIter)
		std::cout << "\nEXIT: Maximum iterations reached.\n";

	// TODO print final stuff

} //generate_trajectory
