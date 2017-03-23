#include "iLQR.h"

#define eps 1e-6 // For finite differencing

// Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
void iLQR::compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u, VecOfMatXd &fx,
												 VecOfMatXd &fu, VecOfVecXd &cx, VecOfVecXd &cu,
												 VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu)
{
    get_dynamics_derivatives(x, u, fx, fu);
    get_cost_derivatives(x, u, cx, cu);
    get_cost_2nd_derivatives(x, u, cxx, cxu, cuu);
}


// -----------------
// Derivative computation

void iLQR::get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
                                        VecOfMatXd &fx, VecOfMatXd &fu)
{
  for (int t=0; t<T; t++)
  {
    fx[t].resize(n,n);
    fu[t].resize(n,m);

    VecXd plus, minus;
    for (int i=0; i<n; i++)
    {
      plus = minus = x[t];
      plus(i) += eps;
      minus(i) -= eps;
      fx[t].col(i) = (integrate_dynamics(plus, u[t])-integrate_dynamics(minus, u[t])) / (2*eps);
    }

    for (int i=0; i<m; i++){
      plus = minus = u[t];
      plus(i) += eps;
      minus(i) -= eps;
      fu[t].col(i) = (integrate_dynamics(x[t], plus)-integrate_dynamics(x[t], minus)) / (2*eps);
    }
    // std::cout << t << '\n';
    // std::cout << "fx[t]: \n" << fx[t] << '\n';
    // std::cout << "fu[t]: \n" << fu[t] << '\n';
    // getchar();
  }
}

void iLQR::get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
															   VecOfVecXd &cx, VecOfVecXd &cu)
{
  for (int t=0; t<T; t++)
  {
    cx[t].resize(n);
    cu[t].resize(m);

    VecXd plus, minus;
    for (int i=0; i<n; i++)
    {
      plus = minus = x[t];
      plus(i) += eps;
      minus(i) -= eps;
      cx[t](i) = (cost(plus, u[t])-cost(minus, u[t])) / (2*eps);
    }

    for (int i=0; i<m; i++){
      plus = minus = u[t];
      plus(i) += eps;
      minus(i) -= eps;
      cu[t](i) = (cost(x[t], plus)-cost(x[t], minus)) / (2*eps);
    }

    // std::cout << t << '\n';
    // std::cout << "cx[t]: \n" << cx[t] << '\n';
    // std::cout << "cu[t]: \n" << cu[t] << '\n';
    // getchar();
  }
}

void iLQR::get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
                              VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu)
{
  for (int t=0; t<T; t++)
  {
    cxx[t].resize(n,n);
    cxu[t].resize(n,m);
    cuu[t].resize(m,m);

    //TODO remove repetition

    VecXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

    //cxx
    for (int i=0; i<n; i++){
      for (int j=i; j<n; j++){
        pp = pm = mp = mm = x[t];
        pp(i) = pm(i) += eps;
        mp(i) = mm(i) -= eps;
        pp(j) = mp(j) += eps;
        pm(j) = mm(j) -= eps;
        cxx[t](i,j) = cxx[t](j,i) = (cost(pp, u[t]) - cost(mp, u[t]) - cost(pm, u[t]) + cost(mm, u[t])) / (4*sqr(eps));
      }
    }
    //cxu
    VecXd px, mx, pu, mu;
    for (int i=0; i<n; i++){
      for (int j=0; j<m; j++){
        px = mx = x[t];
        pu = mu = u[t];
        px(i) += eps;
        mx(i) -= eps;
        pu(j) += eps;
        mu(j) -= eps;
        cxu[t](i,j) = (cost(px, pu) - cost(mx, pu) - cost(px, mu) + cost(mx, mu)) / (4*sqr(eps));
      }
    }
    //cuu
    for (int i=0; i<m; i++){
      for (int j=i; j<m; j++){
        pp = pm = mp = mm = u[t];
        pp(i) = pm(i) += eps;
        mp(i) = mm(i) -= eps;
        pp(j) = mp(j) += eps;
        pm(j) = mm(j) -= eps;
        cuu[t](i,j) = cuu[t](j,i) = (cost(x[t], pp) - cost(x[t], mp) - cost(x[t], pm) + cost(x[t], mm)) / (4*sqr(eps));
      }
    }

    // std::cout << t << '\n';
    // std::cout << "cxx[t]: \n" << cxx[t] << '\n';
    // std::cout << "cxu[t]: \n" << cxu[t] << '\n';
    // std::cout << "cuu[t]: \n" << cuu[t] << '\n';
    // getchar();
  }
} //get_cost_2nd_derivatives
