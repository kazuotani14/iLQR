#include "ilqr.h"
#include "finite_diff.h"

#define eps 1e-6 // For finite differencing

using namespace std::placeholders;

// Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
void iLQR::compute_derivatives(const VecOfVecXd &x, const VecOfVecXd &u)
{
    get_dynamics_derivatives(x, u);
    get_cost_derivatives(x, u);
    get_cost_2nd_derivatives(x, u);
}


// -----------------
// Derivative computation

// Updates fx, fu
void iLQR::get_dynamics_derivatives(const VecOfVecXd &x, const VecOfVecXd &u)
{
  for (int t=0; t<T; t++)
  {
    //TODO figure out how to avoid copying
    VectorXd xi = x[t];
    VectorXd ui = u[t];

    std::function<VectorXd(VectorXd)> dyn_x = std::bind(&Model::dynamics, model, _1, u[t]);
    std::function<VectorXd(VectorXd)> dyn_u = std::bind(&Model::dynamics, model, x[t], _1);

    fx[t] = finite_diff_vec2vec(dyn_x, x[t], 4);
    fu[t] = finite_diff_vec2vec(dyn_u, u[t], 4);
  }
}

// Updates cx, cu
void iLQR::get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u)
{
  for (int t=0; t<T; t++)
  {
    std::function<double(VectorXd)> cost_x = std::bind(&Model::cost, model, _1, u[t]);
    std::function<double(VectorXd)> cost_u = std::bind(&Model::cost, model, x[t], _1);

    cx[t] = finite_diff_vec2scalar(cost_x, x[t]);
    cu[t] = finite_diff_vec2scalar(cost_u, u[t]);
  }
}

// TODO figure out how to reduce repetitive code...
// Update cxx, cxu, cuu
void iLQR::get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u)
{
  int n = x[0].size();
  int m = u[0].size();

  VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  for (int t=0; t<T; t++)
  {
    cxx[t].resize(n,n);
    cxu[t].resize(n,m);
    cuu[t].resize(m,m);

    //cxx
    for (int i=0; i<n; i++){
      for (int j=i; j<n; j++){
        pp = pm = mp = mm = x[t];
        pp(i) = pm(i) += eps;
        mp(i) = mm(i) -= eps;
        pp(j) = mp(j) += eps;
        pm(j) = mm(j) -= eps;
        cxx[t](i,j) = cxx[t](j,i) = (c(pp, u[t]) - c(mp, u[t]) - c(pm, u[t]) + c(mm, u[t])) / (4*sqr(eps));
      }
    }
    //cxu
    VectorXd px, mx, pu, mu;
    for (int i=0; i<n; i++){
      for (int j=0; j<m; j++){
        px = mx = x[t];
        pu = mu = u[t];
        px(i) += eps;
        mx(i) -= eps;
        pu(j) += eps;
        mu(j) -= eps;
        cxu[t](i,j) = (c(px, pu) - c(mx, pu) - c(px, mu) + c(mx, mu)) / (4*sqr(eps));
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
        cuu[t](i,j) = cuu[t](j,i) = (c(x[t], pp) - c(x[t], mp) - c(x[t], pm) + c(x[t], mm)) / (4*sqr(eps));
      }
    }
  }
}
