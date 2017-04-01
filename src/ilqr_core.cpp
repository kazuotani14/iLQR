#include "ilqr.h"
#include "double_integrator.h"
#include "eigen_helpers.h"

#define VERBOSE

double iLQR::init_traj(VecXd &x0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.resize(T+1);
	us.resize(T);

	// call forward_pass, get xs, us, cost
	double cost_i;
	forward_pass(x0, u0, xs, us, cost_i);
	std::cout << "Initial cost: " << cost_i << "\n";
	cost_s = cost_i;

	return cost_s;
}

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
//
void iLQR::generate_trajectory()
{
	// Check initialization
	if (us.size()==0 || xs.size()==0){
		std::cout << "Call init_traj first.\n";
		return;
	}

	VecOfVecXd x = xs;	//nxT
	VecOfVecXd u = us; //2xT

	// Initialize all vectors, matrices we'll be using
	MatrixXd du(2,T); 	//2*T double
	VecOfMatXd fx(T+1); //nxnx(T+1)
	VecOfMatXd fu(T+1); //nx2x(T+1)
	VecOfVecXd cx(T+1); //nx(T+1)
	VecOfVecXd cu(T+1); //2x(T+1)
	VecOfMatXd cxx(T+1); //nxnx(T+1)
	VecOfMatXd cxu(T+1); //nx2x(T+1)
	VecOfMatXd cuu(T+1); //2x2x(T+1)

	Vector2d dV; //2x1
	VecOfVecXd Vx(T+1); //nx(T+1)
	VecOfMatXd Vxx(T+1); //nxnx(T+1)
	VecOfMatXd L(T); //2xnxT
	VecOfVecXd l(T); //2xT

	int n = model->x_dims;
	int m = model->u_dims;
	for (int i=0; i<T+1; i++)
	{
		fx[i] = MatrixXd(n,n);
		fu[i].resize(n,m);
		cx[i].resize(n);
		cu[i].resize(m);
		cxx[i].resize(n,n);
		cxu[i].resize(n,m);
		cuu[i].resize(m,m);
		Vx[i].resize(n);
		Vxx[i].resize(n,n);
		if(i<T)
		{
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
			compute_derivatives(xs,us, fx,fu,cx,cu,cxx,cxu,cuu);
			flgChange = 0;
		}

	} // end top-level for-loop

	std::cout << "end of generate_trajectory" << std::endl;
}

//

// 		// std::cout << "Finished step 1 : compute derivatives. \n";
// 		double t_1 = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
//
// 		//--------------------------------------------------------------------------
// 		// STEP 2: Backward pass, compute optimal control law and cost-to-go
// 		start = std::clock();
//
// 		bool backPassDone = false;
// 		while (!backPassDone)
// 		{
// 	 		// update Vx, Vxx, l, L, dV with back_pass
// 			diverge = 0;
// 			diverge = backward_pass(cx,cu,cxx,cxu,cuu,fx,fu,u, Vx,Vxx,l,L,dV);
// 			// std::cout << "l: \n" << l[0] << '\n';
// 			// std::cout << "L: \n" << L[0] << '\n';
//
// 			if (diverge!=0)
// 			{
// 				// std::cout << "Cholesky failed at timestep " << diverge << ".\n";
// 				dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
// 				lambda    = std::max(lambda * dlambda, lambdaMin);
// 				if (lambda > lambdaMax)
// 						break;
// 				continue;
// 			}
// 			backPassDone = true;
// 		}
// 		// TODO check for termination due to small gradient
// 		double gnorm = 0;
// 		// double gnorm = get_gradient_norm(l, u);
//
// 		// std::cout << "Finished step 2 : backward pass. \n";
//
// 		double t_2 = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
//
// 		//--------------------------------------------------------------------------
// 		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost
// 		start = std::clock();
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
// 	  double t_3 = (std::clock() - start) / (double)(CLOCKS_PER_SEC);
//
// 		std::cout << "Step times: " << t_1 << ' ' << t_2 << ' ' << t_3 << '\n';
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
// } //generate_trajectory
