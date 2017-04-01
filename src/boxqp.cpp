#include "boxqp.h"

/*
	Minimize 0.5*x'*H*x + x'*g  s.t. lower<=x<=upper

	 inputs:
	    H            - positive definite matrix  (m * m)
	    g            - bias vector         			 (m)
	    x0           - initial state             (m)

			lower        - lower bounds              (m)
	    upper        - upper bounds              (m)

	 outputs:
	 		result       - result type (roughly, higher is better, see below)
	    x            - solution = k_i             (m)
	    res.H_free        - subspace cholesky factor=R (n_free * n_free)
	    res.v_free         - set of free dimensions     (m)
											- vector of 0 or 1
*/

// TODO rename (H,g) to (Q,c)

double quadCost(const MatrixXd& H, const VectorXd& g, const VectorXd& x)
{
	return 0.5*x.transpose()*H*x + x.dot(g);
}

VectorXd clamp_to_limits(const VectorXd &x, const VectorXd& lower, const VectorXd& upper)
{
	VectorXd x_clamped(x.size());
  for(int i=0; i<x.size(); i++)
  {
    x_clamped(i) = std::min(upper(i), std::max(lower(i), x(i)));
  }
	return x_clamped;
}

// TODO why not Newton method here?
// Armijo line search: for quadratic cost function with limits on x
// Find a step size in the given search direction that leads to at least the expected decrease in value
lineSearchResult quadclamp_line_search(const VectorXd x0, const VectorXd search_dir,
														 const MatrixXd H, const VectorXd g,
														 const VectorXd lower, const VectorXd upper)
{
	double step = 1;
	lineSearchResult res(x0.size());

	VectorXd grad = H*x0 + g;
	double local_slope = search_dir.dot(grad);
	if(local_slope >= 0) // check if search direction isn't descent direction - this shoudn't happen
	{
		res.failed = true;
		return res;
	}

	VectorXd x_reach = x0 + step*search_dir;
	VectorXd x_clamped = clamp_to_limits(x_reach, lower, upper);
	double v = quadCost(H, g, x_clamped);
	double old_v = quadCost(H, g, x0);

	while ((v - old_v) > Armijo*(step*local_slope))
	{
		step *= stepDec;
		res.n_steps++;

		x_reach = x0 + step*search_dir;
		x_clamped = clamp_to_limits(x_reach, lower, upper);
		v = quadCost(H, g, x_clamped);

		if (step < minStep)
		{
			res.failed = true;
			break;
		}
	}

	res.x_opt = x_clamped;
	res.v_opt = v;
	return res;
}


boxQPResult boxQP(const MatrixXd &H, const VectorXd &g, const VectorXd &x0,
                  const VectorXd& lower, const VectorXd& upper)
{
  int n_dims = x0.size();
  assert(H.cols() == n_dims);
  assert(H.rows() == n_dims);
  assert(g.size() == n_dims);
  assert(lower.size() == n_dims);
  assert(upper.size() == n_dims);

  VectorXd x = clamp_to_limits(x0, lower, upper);
  double val = 0.5*x.transpose()*H*x + x.dot(g);
  double oldvalue = 0;

	int nfactors = 0;
	boxQPResult res(n_dims);

  #ifdef VERBOSE
		std::cout << "==========\nStarting box-QP, dimension " << n_dims << ", initial value: " << val << ".\n";
	#endif

  VectorXd clamped_dims(n_dims);
  VectorXd old_clamped_dims(n_dims);

  for(int iter=0; iter<=qp_maxIter; iter++)
  {
  	if(res.result != 0) break;

    // Check if we've stopped improving
    if(iter>0 && (oldvalue - val) < minRelImprove*std::abs(oldvalue))
    {
      // res.result = 4;
      break;
    }

    VectorXd grad = H*x + g;

    // Find clamped dimensions
    old_clamped_dims = clamped_dims;
    clamped_dims.setZero();
    res.v_free.setOnes();
	std::cout << "checking for clamped" << std::endl;
    for (int i=0; i<n_dims; i++)
    {
	  std::cout << x(i) << " " << grad(i) << std::endl;
      if(x(i)==lower(i) && grad(i)>0)
      {
		std::cout << "clamped" << std::endl;
        clamped_dims(i) = 1;
        res.v_free(i) = 0;
      }
      else if(x(i)==upper(i) && grad(i)<0)
      {
		std::cout << "clamped" << std::endl;
        clamped_dims(i) = 1;
        res.v_free(i) = 0;
      }
    }

    // Check if all dimensions are clamped
    if(clamped_dims.all() == x.size()){
      res.result = 6;
      break;
    }

    // Factorize if clamped dimensions have changed
    if (iter==1 || (old_clamped_dims-clamped_dims).sum() != 0)
    {
      int n_free = res.v_free.sum();
      MatrixXd Hf;

      // TODO remove this hard-coded check - this assumes inputs to boxQP will always be size 2!!
      if (res.v_free[0] == 1)
      {
        Hf = H.block(0, 0, n_free, n_free);
      }
      else
      {
        Hf = H.block(1, 1, n_free, n_free);
      }
      Eigen::LLT<MatrixXd> lltOfHf(Hf); // Cholesky decomposition
      res.H_free = lltOfHf.matrixL().transpose();

      nfactors++;
    }

    // Check if gradient norm is below threshold
    double grad_norm = grad.cwiseProduct(res.v_free).norm();
    if (grad_norm < minGrad){
      // res.result = 5;
      break;
    }

    // get new search direction
    VectorXd grad_clamped = H*(x.cwiseProduct(clamped_dims)) + g;

    VectorXd search = VectorXd::Zero(x.size());

		// TODO remove this hack - assumes size(x)==2
		if(res.v_free[0]==1 && res.v_free[1]==1){
			search = -res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free);
		}
		else if (res.v_free[0]==1){
			search(0) = (-res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free))(0);
		}
		else if (res.v_free[1]==1){
			search(1) = (-res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free))(0);
		}

    lineSearchResult linesearch_res = quadclamp_line_search(x, search, H, g, lower, upper);
    if(linesearch_res.failed)
    {
        res.result = 2;
        break;
    }

    #ifdef VERBOSE
      printf("iter %-3d  value % -9.5g |g| %-9.3g  reduction %-9.3g  linesearch %g^%-2d  n_clamped %d\n",
        iter, linesearch_res.v_opt, grad_norm, oldvalue-linesearch_res.v_opt, stepDec, linesearch_res.n_steps, int(clamped_dims.sum()));
    #endif

    // accept candidate
	x = linesearch_res.x_opt;
	val = linesearch_res.v_opt;

  }

  res.x_opt = x;
  return res;

} //boxQP
