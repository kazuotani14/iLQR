#include "ilqr.h"
#include "double_integrator.h"

#define VERBOSE

double iLQR::init_traj(VectorXd &x0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.resize(T+1);
	us.resize(T);

	// call forward_pass, get xs, us, cost
	double cost_i = forward_pass(x0, u0);
	std::cout << "Initial cost: " << cost_i << "\n";
	cost_s = cost_i;

	//allocate space for later
	du = MatrixXd(2,T);
	fx.resize(T+1);
	fu.resize(T+1);
	cx.resize(T+1);
	cu.resize(T+1);
	cxu.resize(T+1);
	cxx.resize(T+1);
	cuu.resize(T+1);

 	dV = Vector2d(2,1);
	Vx.resize(T+1);
	Vxx.resize(T+1);
	L.resize(T);
	l.resize(T);

	int n = model->x_dims;
	int m = model->u_dims;

	std::fill(fx.begin(), fx.end(), MatrixXd::Zero(n,n));
	std::fill(fu.begin(), fu.end(), MatrixXd::Zero(n,m));
	std::fill(cx.begin(), cx.end(), VectorXd::Zero(n));
	std::fill(cu.begin(), cu.end(), VectorXd::Zero(m));
	std::fill(cxx.begin(), cxx.end(), MatrixXd::Zero(n,n));
	std::fill(cxu.begin(), cxu.end(), MatrixXd::Zero(n,m));
	std::fill(cuu.begin(), cuu.end(), MatrixXd::Zero(m,m));
	std::fill(Vx.begin(), Vx.end(), VectorXd::Zero(n));
	std::fill(Vxx.begin(), Vxx.end(), MatrixXd::Zero(n,n));
	std::fill(L.begin(), L.end(), MatrixXd::Zero(2,n));
	std::fill(l.begin(), l.end(), VectorXd::Zero(2));

	return cost_s;
}

void iLQR::generate_trajectory()
{
	// Check initialization - TODO this shouldn't even be a question... we can just put init_traj here
	if (us.size()==0 || xs.size()==0){
		std::cout << "Call init_traj first.\n";
		return;
	}

	VecOfVecXd x = xs;	//nxT
	VecOfVecXd u = us; //2xT

	// constants, timers, counters
	bool flgChange = true;
	bool stop = false;
	double dcost = 0;
	double z = 0;
	double expected = 0;
	int diverge = 0;

	#ifdef VERBOSE
		std::cout << "\n=========== begin iLQG ===========\n";
	#endif

	int iter;
	for (iter=0; iter<maxIter; iter++)
	{
		if (stop)
			break;

		#ifdef VERBOSE
			std::cout << "Iteration " << iter << ".\n";
		#endif

		//--------------------------------------------------------------------------
		//STEP 1: Differentiate dynamics and cost along new trajectory
		if (flgChange){
			compute_derivatives(xs,us);
			flgChange = 0;
		}
		#ifdef VERBOSE
			std::cout << "Finished step 1 : compute derivatives. \n";
		#endif

		//--------------------------------------------------------------------------
		// STEP 2: Backward pass, compute optimal control law and cost-to-go

		bool backPassDone = false;
		while (!backPassDone)
		{
	 		// update Vx, Vxx, l, L, dV with back_pass
			diverge = 0;
			diverge = backward_pass();

			if (diverge != 0)
			{
				#ifdef VERBOSE
					std::cout << "Cholesky failed at timestep " << diverge << ".\n";
				#endif

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

		#ifdef VERBOSE
			std::cout << "Finished step 2 : backward pass. \n";
		#endif


	} // end top-level for-loop

	std::cout << "end of generate_trajectory" << std::endl;

//

//

////
// 		//--------------------------------------------------------------------------
// 		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost
//
// 		bool fwdPassDone = 0;
// 		VecOfVecXd xnew(trajectoryLength+1);
// 		VecOfVecXd unew(trajectoryLength);
// 		double new_cost, alpha;
//
// 		if (backPassDone) //  serial backtracking line-search
// 		{
// 			for (int i=0; i<Alpha.size(); i++){
// 				alpha = Alpha(i);
// 				forward_pass(x_0, adjust_u(us,l,alpha), xnew,unew,new_cost, x,L);
// 				dcost    = cost_s - new_cost;
// 				expected = -alpha * (dV(0) + alpha*dV(1));
//
// 				// std::cout << "dV: " << dV(0) << ' ' << dV(1) << '\n';
// 				// std::cout << "dcost: " << dcost << '\n';
// 				// std::cout << "expected: " << expected << '\n';
// 				// getchar();
//
// 				if (expected>0){
// 					z = dcost/expected;
// 				}
// 				else{
// 					z = sgn(dcost);
// 					std::cout << "Warning: non-positive expected reduction: should not occur\n";
// 				}
// 				if(z > zMin){
// 					fwdPassDone = 1;
// 					break;
// 				}
// 			}
//
// 			if (!fwdPassDone){
// 				alpha = 0.0; // signals failure of forward pass
// 			}
// 		}
//
// 		// std::cout << "Finished step 3 : forward pass. \n";
//
// 	//--------------------------------------------------------------------------
// 	// STEP 4: accept step (or not), print status
// 		#ifdef VERBOSE
// 	 	if (iter==0)
// 	 		std::cout << "iteration\tcost\t\treduction\texpected\tgradient\tlog10(lambda)\n";
// 		#endif
//
// 	 	if (fwdPassDone)
// 	 	{
// 			#ifdef VERBOSE
// 				printf("%-12d\t%-12.6g\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
// 	                iter, new_cost, dcost, expected, gnorm, log10(lambda));
// 			#endif
//
// 	 		// decrease lambda
// 	 		dlambda   = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
// 	 		lambda    = lambda * dlambda * (lambda > lambdaMin);
//
// 	 		// accept changes
// 	 		us              = unew;
// 	 		xs              = xnew;
// 	 		cost_s          = new_cost;
// 	 		flgChange       = true;
//
// 	 		//terminate?
// 	 		if (dcost < tolFun){
// 				#ifdef VERBOSE
// 		 			std::cout << "\nSUCCESS: cost change < tolFun\n";
// 				#endif
// 	 			break;
// 	 		}
// 	 	}
// 	 	else // no cost improvement
// 	 	{
// 	 		// increase lambda
// 	 		dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
// 	 		lambda   = std::max(lambda * dlambda, lambdaMin);
//
// 			#ifdef VERBOSE
// 				printf("%-12d\t%-12s\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
// 	                iter,"NO STEP", dcost, expected, gnorm, log10(lambda));
// 			#endif
//
// 	 		// terminate?
// 	 		if (lambda > lambdaMax){
// 				#ifdef VERBOSE
// 	 				std::cout << "\nEXIT: lambda > lambdaMax\n";
// 				#endif
// 	 			break;
// 	 		}
// 		}
// 	} // optimization for loop
//
// 	#ifdef VERBOSE
// 		if (iter==maxIter) std::cout << "\nEXIT: Maximum iterations reached.\n";
// 	#endif
//
// 	// output_to_csv();
//
} //generate_trajectory


// VecOfVecXd iLQR::adjust_u(VecOfVecXd &u, VecOfVecXd &l, double alpha)
// {
// 	VecOfVecXd new_u = u;
// 	for(int i=0; i<u.size(); i++){
// 		new_u[i] += l[i]*alpha;
// 	}
// 	return new_u;
// }
//
// void iLQR::output_to_csv()
// {
// 	FILE *XU = fopen("XU.csv", "w");
// 	for(int t=0; t<T; t++) {
// 			fprintf(XU, "%f, %f, %f, %f, %f, %f, ",
// 									xs[t](0), xs[t](1), xs[t](2), xs[t](3), xs[t](4), xs[t](5));
// 			fprintf(XU, "%f, %f \n", us[t](0), us[t](1));
// 	}
// 	fclose(XU);
// }
//
// // double iLQR::get_gradient_norm(VecOfVecXd l, VecOfVecXd u)
// // {
// // 	for (int i=0; i<l.size())
// // }
