#include "iLQR.h"

#define VERBOSE 0

VecXd iLQR::clamp_to_limits(VecXd &u)
{
	VecXd u_clamped(m);
	u_clamped(0) = std::min(control_limits(0,1), std::max(control_limits(0,0), u(0)));
	u_clamped(1) = std::min(control_limits(1,1), std::max(control_limits(1,0), u(1)));
	return u_clamped;
}

VecXd iLQR::subvec_w_ind(VecXd &vec, VecXd &indices)
{
	VecXd subvec;
	for(int i=0; i<vec.size(); i++){
		if(indices(i)>0){
			push_back(subvec, vec(i));
		}
	}
	return subvec;
}

int iLQR::boxQP(MatXd &H, VecXd &g, VecXd &x0,  VecXd &x, MatXd& Hfree, VecXd &free_v)
{
/*
	Minimize 0.5*x'*H*x + x'*g  s.t. lower<=x<=upper

	 inputs:
	    H            - positive definite matrix from QuuF		 (m * m)
	    g            - bias vector from Qu        					 (m)
	    x0           - initial state              (m)

			lower        - lower bounds from control_limits      (m)
	    upper        - upper bounds from control_limit       (m)

	 outputs:
	 		=result       - result type (roughly, higher is better, see below)
	    x            - solution = k_i             (m)
	    Hfree        - subspace cholesky factor=R (n_free * n_free)
	    free_v         - set of free dimensions     (m)
											- vector of 0 or 1
*/

	free_v.resize(m);

	// Optimization parameters
	int maxIter        		= 100;      // maximum number of iterations
	double minGrad        = 1e-8;     // minimum norm of non-fixed gradient
	double minRelImprove  = 1e-8;     // minimum relative improvement
	double stepDec        = 0.6;     	// factor for decreasing stepsize
	double minStep        = 1e-22;    // minimal stepsize for linesearch
	double Armijo         = 0.1;   		// Armijo parameter (fraction of linear improvement required)

	VecXd lower = control_limits.col(0);
	VecXd upper = control_limits.col(1);

	// initial state
	x = clamp_to_limits(x0);

	// initial objective value
	double value = x.dot(g) + 0.5*x.transpose()*H*x;

	if (VERBOSE){
		std::cout << "==========\nStarting box-QP, dimension " << m << ", initial value: " << value << ".\n";
	}

	int result = 0;
	int nfactor = 0;
	double oldvalue = 0;
	bool factorize = false;

	VecXd clamped(m);
	VecXd old_clamped(m);
	for (int iter=1; iter<=maxIter; iter++)
	{
		if (result != 0)
			break;

		// check relative improvement
		if (iter>1 && (oldvalue - value)<minRelImprove*std::abs(oldvalue)) {
			result = 4;
			break;
		}
		oldvalue = value;

		// Get gradient
		VecXd grad = g + H*x;

		// Find clamped dimensions
		old_clamped = clamped;
		clamped.setZero();
		free_v.setOnes();
		for (int i=0; i<m; i++){
			if(x(i)==lower(i) && grad(i)>0){
				clamped(i) = 1;
				free_v(i) = 0;
			}
			else if(x(i)==upper(i) && grad(i)<0){
				clamped(i) = 1;
				free_v(i) = 0;
			}
		}

		// Check for all clamped TODO
		if(clamped.all() == m){
			result = 6;
			break;
		}

		// Factorize if clamped has changed
		if (iter==1){
			factorize = true;
		}
		else if ((old_clamped-clamped).sum() != 0){
			factorize = true;
		}
		else{
			factorize = false;
		}

		if (factorize){
			int n_free = free_v.sum();
			MatXd Hf;
			// TODO fix this hard-coded check
			if (free_v[0] == 1){
				Hf = H.block(0, 0, n_free, n_free);
			}
			else{
				Hf = H.block(1, 1, n_free, n_free);
			}
			Eigen::LLT<MatXd> lltOfHf(Hf); // compute the Cholesky decomposition of A
			Hfree = lltOfHf.matrixL().transpose();

			nfactor++;
		}

		// check gradient norm
		double gnorm = grad.cwiseProduct(free_v).norm();
		if (gnorm < minGrad){
			result = 5;
			break;
		}

		// get search direction
		VecXd grad_clamped = g + H*(x.cwiseProduct(clamped));

		VecXd search(m);
		search.setZero();
		// TODO fix this hack
		if(free_v[0]==1 && free_v[1]==1){
			search = -Hfree.inverse() * (Hfree.transpose().inverse()*subvec_w_ind(grad_clamped, free_v)) - subvec_w_ind(x, free_v);
		}
		else if (free_v[0]==1){
			search(0) = (-Hfree.inverse() * (Hfree.transpose().inverse()*subvec_w_ind(grad_clamped, free_v)) - subvec_w_ind(x, free_v))(0);
		}
		else if (free_v[1]==1){
			search(1) = (-Hfree.inverse() * (Hfree.transpose().inverse()*subvec_w_ind(grad_clamped, free_v)) - subvec_w_ind(x, free_v))(0);
		}

		// check for descent direction
		double sdotg = (search.cwiseProduct(grad)).sum();
		if (sdotg >= 0)
			break; // This shouldn't happen


		// Armijo line search
		double step = 1;
		int nstep = 0;
		VecXd reach = x + step*search;
		VecXd xc = clamp_to_limits(reach);
		double vc = xc.dot(g) + 0.5*xc.dot(H*xc); // TODO use dot?

		while ((vc - oldvalue)/(step*sdotg) < Armijo){
			step *= stepDec;
			nstep++;
			reach = x + step*search;
			xc = clamp_to_limits(reach);
			vc = xc.dot(g) + 0.5*xc.dot(H*xc);
			if (step<minStep){
				result = 2;
				break;
			}
		}

		if (VERBOSE){
			printf("iter %-3d  value % -9.5g |g| %-9.3g  reduction %-9.3g  linesearch %g^%-2d  n_clamped %d\n",
				iter, vc, gnorm, oldvalue-vc, stepDec, nstep, int(clamped.sum()));
		}

		// accept candidate
		x = xc;
		value = vc;

} //iterations

	return result;

} //boxQP

void iLQR::demoQP()
{
	// For debugging
	VecXd g(2);
	g << -10, 99;
	MatXd H(2,2);
	H << 0.1, 0.2, 0.3, -44;
	H = H*H.transpose();
	Vec2d lower_new(-1,-1);
	Vec2d upper_new(1,1);
	control_limits.col(0) = lower_new;
	control_limits.col(1) = upper_new;

	VecXd x0(2);
	x0 << -1.5, 0;

	VecXd x;
	MatXd Hfree;
	VecXd free_v;

	int result = boxQP(H, g, x0,  x, Hfree, free_v);
	std::cout << "\n--------\n";
	std::cout << "x: "; print_vec(x);
	std::cout << "result: " << result << std::endl;
	std::cout << "Hfree: " << Hfree << std::endl;
	std::cout <<  "free: "; print_vec(free_v);
}
