#include "iLQR2.h"

#define eps 0.0001

// This has a weird condition based on number of outputs in MATLAB code
// Split into two functions
double iLQR::get_dynamics_and_cost(const VecXd &x, const VecXd u, VecXd &x1)
{
	// returns cost, modifies dx
		x1 = dynamics(x,u);
		double c = cost(x,u);
		return c;
}

void iLQR::compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd &fx,
												 VecOfMatXd &fu, VecOfVecXd &cx, VecOfVecXd &cu,
												 VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu)
{
  for (int t=0; t<T; t++)
  {
    get_dynamics_derivatives(x[t], u[t], fx[t], fu[t]);
    get_cost_derivatives(x[t], u[t], cx[t], cu[t]);
    get_cost_2nd_derivatives(x[t], u[t], cxx[t], cxu[t], cuu[t]);
  }
}


// -----------------
// Derivative computation

void iLQR::get_dynamics_derivatives(const VecXd &x, const VecXd &u,
															MatXd &fx, MatXd &fu)
{
  VecXd plus, minus;
  for (int i=0; i<n; i++)
  {
    plus = minus = x;
    plus(i) += eps;
    minus(i) += eps;
    fx.col(i) = (dynamics(plus, u)-dynamics(minus, u)) / (2*eps);
  }

  for (int i=0; i<m; i++){
    plus = minus = u;
    plus(i) += eps;
    minus(i) += eps;
    fu.col(i) = (dynamics(x, plus)-dynamics(x, minus)) / (2*eps);
  }
}

void iLQR::get_cost_derivatives(const VecXd &x, const VecXd &u,
															VecXd &cx, VecXd &cu)
{
  VecXd plus, minus;
  for (int i=0; i<n; i++)
  {
    plus = minus = x;
    plus(i) += eps;
    minus(i) += eps;
    cx(i) = (cost(plus, u)-cost(minus, u)) / (2*eps);
  }

  for (int i=0; i<m; i++){
    plus = minus = u;
    plus(i) += eps;
    minus(i) += eps;
    cu(i) = (cost(x, plus)-cost(x, minus)) / (2*eps);
  }
}

void iLQR::get_cost_2nd_derivatives(const VecXd &x, const VecXd &u,
															MatXd &cxx, MatXd &cxu, MatXd &cuu)
{
  //TODO remove repetition
  VecXd pp, pm, mp, mm; //plus-plus, plus-minus, ....
  //cxx
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      pp = pm = mp = mm = x;
      pp(i) += eps; pp(j) += eps;
      pm(i) += eps; pm(j) += eps;
      mp(i) -= eps; mp(j) -= eps;
      mm(i) -= eps; mm(j) -= eps;
      cxx(i,j) = (cost(mm, u) + cost(pp, u) + cost(pm, u) + cost(mp, u)) / (4*sqr(eps));
    }
  }
  //cxu
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++){
      pp = pm = mp = mm = x;
      pp(i) += eps; pp(j) += eps;
      pm(i) += eps; pm(j) += eps;
      mp(i) -= eps; mp(j) -= eps;
      mm(i) -= eps; mm(j) -= eps;
      cxu(i,j) = (cost(mm, u) + cost(pp, u) + cost(pm, u) + cost(mp, u)) / (4*sqr(eps));
    }
  }
  //cuu
  for (int i=0; i<m; i++){
    for (int j=0; j<m; j++){
      pp = pm = mp = mm = x;
      pp(i) += eps; pp(j) += eps;
      pm(i) += eps; pm(j) += eps;
      mp(i) -= eps; mp(j) -= eps;
      mm(i) -= eps; mm(j) -= eps;
      cuu(i,j) = (cost(mm, u) + cost(pp, u) + cost(pm, u) + cost(mp, u)) / (4*sqr(eps));
    }
  }
} //get_cost_2nd_derivatives
