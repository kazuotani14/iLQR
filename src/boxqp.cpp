#include "boxqp.h"

/*
  Minimize 0.5*x'*Q*x + x'*c  s.t. lower<=x<=upper
  See Appendix in [Tassa 2014]

  What's cool about this QP solver vs. generic ones: 
  If initial guess has same active set/clamped dimensions as optimum, it will be reached in one iteration!

   inputs:
      Q            - positive definite matrix  (m * m)
      c            - bias vector                (m)
      x0           - initial guess             (m)

      lower        - lower bounds              (m)
      upper        - upper bounds              (m)

   outputs:
     result       - result type (roughly, higher is better, see below)
     x            - solution = k_i             (m)
     res.R_free   - subspace cholesky factor=R (n_free * n_free)
     res.v_free   - set of free dimensions     (m)
                      - vector of 0 or 1
*/

boxQPResult boxQP(const MatrixXd& Q, const VectorXd& c, const VectorXd& x0,
                  const VectorXd& lower, const VectorXd& upper) {
  int n_dims = x0.size();
  assert(Q.cols() == n_dims);
  assert(Q.rows() == n_dims);
  assert(c.size() == n_dims);
  assert(lower.size() == n_dims);
  assert(upper.size() == n_dims);

  VectorXd x = clamp_to_limits(x0, lower, upper);
  double val = x.transpose()*Q*x + x.dot(c);
  double oldvalue = 0;

  int nfactors = 0;
  boxQPResult res(n_dims);

  #ifdef VERBOSE
    std::cout << "==========\nStarting box-QP, dimension " << n_dims << ", initial value: " << val << ".\n";
  #endif

  VectorXd clamped_dims(n_dims); // TODO change this to ints - need to deal with casting when interacting with other VectorXd's
  VectorXd old_clamped_dims(n_dims);
  VectorXd grad(n_dims), grad_clamped(n_dims), search(n_dims);

  for(int iter=0; iter<=qp_maxIter; iter++) {
    if(res.result != 0) break; // TODO we might not need

    // Check if we've stopped improving
    if(iter>0 && (oldvalue - val) < minRelImprove*std::abs(oldvalue)) {
      res.result = 4;
      break;
    }
    grad = Q*x + c;
    oldvalue = val;

    // Find clamped dimensions
    old_clamped_dims = clamped_dims;
    clamped_dims.setZero();
    res.v_free.setOnes();
    for (int i=0; i<n_dims; i++) {
      // approx_eq instead of equals to deal with numerical issues
      if( (approx_eq(x(i), lower(i)) && grad(i)>0) || (approx_eq(x(i), upper(i)) && grad(i)<0) ) {
        clamped_dims(i) = 1;
        res.v_free(i) = 0;
      }
    }

    // Check if all dimensions are clamped. 
    if(clamped_dims.all()) {
      res.result = 6;
      break;
    }

    // Factorize if clamped dimensions have changed
    if (iter==0 || (old_clamped_dims-clamped_dims).sum() != 0) {
      // Cholesky decomposition - so taking inverse later is more efficient
      MatrixXd Qfree;
      Qfree = extract_bool_rowsandcols(Q, res.v_free);

      Eigen::LLT<MatrixXd> choleskyOfQfree(Qfree); 
      if(choleskyOfQfree.matrixL().size() > 0) { // I'm not sure why this happens, but just in case TODO debug
        res.R_free = choleskyOfQfree.matrixL().transpose();
      }
      nfactors++;
    }

    // Check if gradient norm is below threshold
    double grad_norm = subvec_w_ind(grad, res.v_free).norm();
    if (grad_norm < minGrad) {
      res.result = 5;
      break;
    }

    // get new search direction
    grad_clamped = Q*(x.cwiseProduct(clamped_dims)) + c;

    // search(free) = -Hfree \ (Hfree' \ grad_clamped(free)) - x(free);
    search.setZero();
    if (res.v_free.all()) {
      search = -res.R_free.inverse() * res.R_free.transpose().inverse()
               * subvec_w_ind(grad_clamped, res.v_free)
               - subvec_w_ind(x, res.v_free);
    }
    else {
      VectorXd search_update_vals = -(res.R_free.inverse() * res.R_free.transpose().inverse())
                                    *subvec_w_ind(grad_clamped, res.v_free)
                                    - subvec_w_ind(x, res.v_free);
      
      // TODO find cleaner way
      int update_idx = 0;
      for (int k=0; k<search.size(); ++k) {
        if(res.v_free[k]==1) search(k) = search_update_vals(update_idx++);
      }
    }

    lineSearchResult linesearch_res = quadclamp_line_search(x, search, Q, c, lower, upper);
    if(linesearch_res.failed) {
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

// Armijo line search: for quadratic cost function with limits on x
// Find a step size in the given search direction that leads to at least the expected decrease in value
lineSearchResult quadclamp_line_search(const VectorXd& x0, const VectorXd& search_dir,
                     const MatrixXd& Q, const VectorXd& c,
                     const VectorXd& lower, const VectorXd& upper) {
  double step = 1;
  lineSearchResult res(x0.size());

  VectorXd grad = Q*x0 + c;
  double local_slope = search_dir.dot(grad);
  if(local_slope >= 0) { // check if search direction isn't descent direction - this shoudn't happen
    res.failed = true;
    return res;
  }

  VectorXd x_reach = x0 + step*search_dir;
  VectorXd x_clamped = clamp_to_limits(x_reach, lower, upper);
  double v = quadCost(Q, c, x_clamped);
  double old_v = quadCost(Q, c, x0);

  while ((v - old_v)/(step*local_slope) < Armijo) {
    step *= stepDec;
    res.n_steps++;

    x_reach = x0 + step*search_dir;
    x_clamped = clamp_to_limits(x_reach, lower, upper);
    v = quadCost(Q, c, x_clamped);

    if (step < minStep) {
      res.failed = true;
      break;
    }
  }

  res.x_opt = x_clamped;
  res.v_opt = v;
  return res;
}
