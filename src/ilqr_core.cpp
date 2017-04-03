#include "ilqr.h"
#include "double_integrator.h"

#define SHOWPROGRESS
// #define VERBOSE

double iLQR::init_traj(VectorXd &x_0, VecOfVecXd &u0)
{
	//initialize xs and us, return or store initial cost
	xs.resize(T+1);
	us.resize(T);
	x0 = x_0;

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
	k.resize(T);
	K.resize(T);

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
	std::fill(k.begin(), k.end(), VectorXd::Zero(m));
	std::fill(K.begin(), K.end(), MatrixXd::Zero(m,n));

	return cost_s;
}

void iLQR::generate_trajectory()
{
	// Check initialization - TODO this shouldn't even be a question... we can just put init_traj here
	if (us.size()==0 || xs.size()==0){
		std::cout << "Call init_traj first.\n";
		return;
	}

	// TODO check if we ever use this
	VecOfVecXd x_old, u_old;

	// constants, timers, counters
	bool flgChange = true;
	bool stop = false;
	double dcost = 0;
	double z = 0;
	double expected = 0;
	int diverge = 0;
	double new_cost, gnorm;

	#ifdef VERBOSE
		std::cout << "\n=========== begin iLQG ===========\n";
	#endif

	int iter;
	for (iter=0; iter<maxIter; iter++)
	{
		x_old = xs; u_old = us;

		if (stop)
			break;

		#ifdef VERBOSE
			std::cout << "Iteration " << iter << ".\n";
		#endif

		//--------------------------------------------------------------------------
		//STEP 1: Differentiate dynamics and cost along new trajectory

		auto start = std::chrono::system_clock::now();

		if (flgChange)
		{
			compute_derivatives(xs,us);
			flgChange = 0;
		}
		#ifdef VERBOSE
			cout << "Finished step 1 : compute derivatives." << endl;;
		#endif

		auto now = std::chrono::system_clock::now();
		long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
		cout << "compute_derivatives took: " << elapsed/1000. << " seconds." << endl;

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
					std::cout << "Backpass failed at timestep " << diverge << ".\n";
				#endif

				dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
				lambda    = std::max(lambda * dlambda, lambdaMin);
				if (lambda > lambdaMax)
						break;
				continue;
			}
			backPassDone = true;
		}

		// check for termination due to small gradient
		gnorm = get_gradient_norm(k, us);
		if (gnorm < tolGrad && lambda < 1e-5)
		{
			#ifdef SHOWPROGRESS
				cout << "\nSUCCESS: gradient norm < tolGrad\n" << endl;
			#endif
			break;
		}

		#ifdef VERBOSE
			std::cout << "Finished step 2 : backward pass. \n";
		#endif

		//--------------------------------------------------------------------------
		// STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost

		bool fwdPassDone = 0;
		VecOfVecXd xnew(T+1);
		VecOfVecXd unew(T);
		double alpha;

		if (backPassDone) //  serial backtracking line-search
		{
			for (int i=0; i<Alpha.size(); i++)
			{
				alpha = Alpha(i);
				VecOfVecXd u_plus_feedforward = add_bias_to_u(us, k, alpha);

				new_cost = forward_pass(x0, u_plus_feedforward);
				dcost    = cost_s - new_cost;
				expected = -alpha * (dV(0) + alpha*dV(1));

				if (expected>0)
				{
					z = dcost/expected;
				}
				else
				{
					z = sgn(dcost);
					#ifdef VERBOSE
						cout << "Warning: non-positive expected reduction: should not occur" << endl;
					#endif
				}

				if(z > zMin)
				{
					fwdPassDone = 1;
					break;
				}
			}

			if (!fwdPassDone)
			{
				#ifdef VERBOSE
					cout << "Forward pass failed" << endl;
				#endif
				alpha = 0.0; // signals failure of forward pass
			}
		}

		#ifdef VERBOSE
			cout << "Finished step 3 : forward pass." <<endl;
		#endif

		//--------------------------------------------------------------------------
		// STEP 4: accept step (or not), print status
		#ifdef SHOWPROGRESS
			if (iter==0)
	 		std::cout << "iteration\tcost\t\treduction\texpect\t\tgrad\t\tlog10(lambda)\n";
		#endif

	 	if (fwdPassDone)
	 	{
			#ifdef SHOWPROGRESS
				printf("%-12d\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
	                iter, new_cost, dcost, expected, gnorm, log10(lambda));
			#endif

	 		// decrease lambda
	 		dlambda   = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
	 		lambda    = lambda * dlambda * (lambda > lambdaMin);

	 		// accept changes
			// cout << "accepting new cost: " << new_cost << endl;
	 		cost_s          = new_cost;
	 		flgChange       = true;

	 		//terminate?
	 		if (dcost < tolFun)
			{
				#ifdef SHOWPROGRESS
		 			std::cout << "\nSUCCESS: cost change < tolFun\n";
				#endif
	 			break;
	 		}
	 	}
	 	else // no cost improvement
	 	{
			xs = x_old;
			us = u_old;

	 		// increase lambda
	 		dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
	 		lambda   = std::max(lambda * dlambda, lambdaMin);

			#ifdef SHOWPROGRESS
				printf("%-12d\t%-12s\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
	                iter,"NO STEP", dcost, expected, gnorm, log10(lambda));
			#endif

	 		// terminate?
	 		if (lambda > lambdaMax){
				#ifdef SHOWPROGRESS
	 				std::cout << "\nEXIT: lambda > lambdaMax\n";
				#endif
	 			break;
	 		}
		}

		#ifdef SHOWPROGRESS
			if (iter==maxIter) std::cout << "\nEXIT: Maximum iterations reached.\n";
		#endif

	} // end top-level for-loop

	output_to_csv("ilqr_result.csv");

} //generate_trajectory

VecOfVecXd iLQR::add_bias_to_u(const VecOfVecXd &u, const VecOfVecXd &l, const double alpha)
{
	VecOfVecXd new_u = u;
	for(int i=0; i<u.size(); i++){
		new_u[i] += l[i]*alpha;
	}
	return new_u;
}

// Replaces this line from matlab: g_norm = mean(max(abs(l) ./ (abs(u)+1), [], 1));
double iLQR::get_gradient_norm(const VecOfVecXd& l, const VecOfVecXd& u)
{
	std::vector<double> vals(l.size());
	for (int i=0; i<l.size(); i++)
	{
		VectorXd v = l[i].cwiseAbs().array() / (u[i].cwiseAbs().array() + 1);
		double max_val = v.maxCoeff();
		vals[i] = max_val;
	}
	double average = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
	return average;
}

void iLQR::output_to_csv(const std::string filename)
{
	FILE *XU = fopen(filename.c_str(), "w");

	for(int i=1; i<=xs[0].size(); i++)
		fprintf(XU, "x%d, ", i);
	for(int j=0; j<us[0].size(); j++)
		fprintf(XU, "u%d, ", j);
	fprintf(XU, "u%d\n", int(us[0].size()));

	for(int t=0; t<T; t++)
	{
		for(int i=0; i<xs[t].size(); i++)
			fprintf(XU, "%f, ", xs[t](i));
		for(int j=0; j<us[t].size()-1; j++)
			fprintf(XU, "%f, ", us[t](j));
		fprintf(XU, "%f\n", us[t](us[t].size()-1));
	}

	for(int i=0; i<xs[T].size(); i++)
		fprintf(XU, "%f, ", xs[T](i));

	fclose(XU);
	cout << "Saved iLQR result to " << filename << endl;
}
