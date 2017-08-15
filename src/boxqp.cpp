#include "boxqp.h"

/*
  Minimize 0.5*x'*Q*x + x'*c  s.t. lower<=x<=upper

   inputs:
      Q            - positive definite matrix  (m * m)
      c            - bias vector                (m)
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

double quadCost(const MatrixXd& Q, const VectorXd& c, const VectorXd& x)
{
  return 0.5*x.transpose()*Q*x + x.dot(c);
}

VectorXd clamp_to_limits(const VectorXd &x, const VectorXd& lower, const VectorXd& upper)
{
  VectorXd x_clamped(x.size());
  for(int i=0; i<x.size(); i++)
  {
    x_clamped(i) = std::min(upper(i), std::max(lower(i), x(i)));
  }
  // VectorXd x_clamped = upper.cwiseMin(x.cwiseMax(lower));
  return x_clamped;
}

// Armijo line search: for quadratic cost function with limits on x
// Find a step size in the given search direction that leads to at least the expected decrease in value
lineSearchResult quadclamp_line_search(const VectorXd& x0, const VectorXd& search_dir,
                     const MatrixXd& Q, const VectorXd& c,
                     const VectorXd& lower, const VectorXd& upper)
{
  double step = 1;
  lineSearchResult res(x0.size());

  VectorXd grad = Q*x0 + c;
  double local_slope = search_dir.dot(grad);
  if(local_slope >= 0) // check if search direction isn't descent direction - this shoudn't happen
  {
    res.failed = true;
    return res;
  }

  VectorXd x_reach = x0 + step*search_dir;
  VectorXd x_clamped = clamp_to_limits(x_reach, lower, upper);
  double v = quadCost(Q, c, x_clamped);
  double old_v = quadCost(Q, c, x0);

  while ((v - old_v)/(step*local_slope) < Armijo)
  {
    step *= stepDec;
    res.n_steps++;

    x_reach = x0 + step*search_dir;
    x_clamped = clamp_to_limits(x_reach, lower, upper);
    v = quadCost(Q, c, x_clamped);

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


boxQPResult boxQP(const MatrixXd &Q, const VectorXd &c, const VectorXd &x0,
                  const VectorXd& lower, const VectorXd& upper)
{
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

  VectorXd clamped_dims(n_dims);
  VectorXd old_clamped_dims(n_dims);

  for(int iter=0; iter<=qp_maxIter; iter++)
  {
    if(res.result != 0) break;

    // Check if we've stopped improving
    if(iter>0 && (oldvalue - val) < minRelImprove*std::abs(oldvalue))
    {
      res.result = 4;
      break;
    }
  VectorXd grad = Q*x + c;
  oldvalue = val;

    // Find clamped dimensions
    old_clamped_dims = clamped_dims;
    clamped_dims.setZero();
    res.v_free.setOnes();
    for (int i=0; i<n_dims; i++)
    {
      if(approx_eq(x(i), lower(i)) && grad(i)>0)
      {
        clamped_dims(i) = 1;
        res.v_free(i) = 0;
      }
      else if(approx_eq(x(i), upper(i)) && grad(i)<0)
      {
        clamped_dims(i) = 1;
        res.v_free(i) = 0;
      }
    }

    // Check if all dimensions are clamped
    if(clamped_dims.all())
  {
      res.result = 6;
      break;
    }

    // Factorize if clamped dimensions have changed
    if (iter==0 || (old_clamped_dims-clamped_dims).sum() != 0)
    {
      int n_free = res.v_free.sum();
      MatrixXd Qfree;

      // TODO remove this hard-coded check - this assumes inputs to boxQP will always be size 2!!
    // what we want is: Qfree = Q(v_free, v_  free)
      if (res.v_free[0] == 1)
      {
        Qfree = Q.block(0, 0, n_free, n_free);
      }
      else
      {
        Qfree = Q.block(1, 1, n_free, n_free);
      }

      Eigen::LLT<MatrixXd> choleskyOfQfree(Qfree); // Cholesky decomposition
    if(choleskyOfQfree.matrixL().size() > 0) // I don't know why this happens...
        res.H_free = choleskyOfQfree.matrixL().transpose();

      nfactors++;
    }

    // Check if gradient norm is below threshold
    double grad_norm = grad.cwiseProduct(res.v_free).norm();
    if (grad_norm < minGrad)
  {
      res.result = 5;
      break;
    }

    // get new search direction
    VectorXd grad_clamped = Q*(x.cwiseProduct(clamped_dims)) + c;

    VectorXd search = VectorXd::Zero(x.size());

  // TODO remove this hack - assumes size(x)==2
  if(res.v_free[0]==1 && res.v_free[1]==1)
  {
    search = -res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free);
  }
  else if (res.v_free[0]==1){
    search(0) = (-res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free))(0);
  }
  else if (res.v_free[1]==1){
    search(1) = (-res.H_free.inverse() * (res.H_free.transpose().inverse()*subvec_w_ind(grad_clamped, res.v_free)) - subvec_w_ind(x, res.v_free))(0);
  }

  lineSearchResult linesearch_res = quadclamp_line_search(x, search, Q, c, lower, upper);
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
