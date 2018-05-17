#include "ilqr.h"
#include "finite_diff.h"
#include <thread>
#include <functional>

#include "omp.h"
#define EIGEN_DONT_PARALLELIZE
#define OMP_NUM_THREADS 4

#define eps2 1e-3 // For finite differencing

using namespace std::placeholders;

// Updates fx, fu
void iLQR::get_dynamics_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& f_x, VecOfMatXd& f_u) {
  int state_size = x[0].size();

  //#pragma omp parallel for 
  for (int t=0; t<T; t++) {
    std::function<VectorXd(VectorXd)> dyn_x = std::bind(&Model::integrate_dynamics, model, _1, u[t], dt);
    std::function<VectorXd(VectorXd)> dyn_u = std::bind(&Model::integrate_dynamics, model, x[t], _1, dt);

    f_x[t] = finite_diff_jacobian(dyn_x, x[t], state_size); 
    f_u[t] = finite_diff_jacobian(dyn_u, u[t], state_size);
  }
}

// Updates cx, cu
void iLQR::get_cost_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfVecXd& c_x, VecOfVecXd& c_u) {
  int u_dims = u[0].size();

  //#pragma omp parallel for
  for (int t=0; t<T+1; t++) {
    VectorXd ut(u_dims);
    if(t<T) 
      ut = u[t];
    else 
      ut.setZero();

    std::function<double(VectorXd)> cost_x = std::bind(&Model::cost, model, _1, ut);
    std::function<double(VectorXd)> cost_u = std::bind(&Model::cost, model, x[t], _1);
    std::function<double(VectorXd)> cost_f = std::bind(&Model::final_cost, model, _1);

    if(t<T) {
      c_x[t] = finite_diff_gradient(cost_x, x[t]);
      c_u[t] = finite_diff_gradient(cost_u, ut);
    }
    else {
      c_x[t] = finite_diff_gradient(cost_f, x[t]);
      c_u[t].resize(u_dims);
      c_u[t].setZero();
    }
  }
}

// Update cxx, cxu, cuu
void iLQR::get_cost_2nd_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xx, VecOfMatXd& c_xu, VecOfMatXd& c_uu)
{
  VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  int x_dims = x[0].size();
  int u_dims = u[0].size();
  for (int t=0; t<T+1; t++) {
    c_xx[t].resize(x_dims, x_dims);
    c_xu[t].resize(x_dims, u_dims);
    c_uu[t].resize(u_dims, u_dims); // TODO why is this calculated up to T? shouldn't it end at T-1?
  } 

  calculate_cxx(x, u, c_xx, 0, T+1);
  calculate_cxu(x, u, c_xu, 0, T+1);
  calculate_cuu(x, u, c_uu, 0, T+1);
}

void iLQR::calculate_cxx(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xx, int start_T, int end_T) {
  std::function<double(VectorXd)> c;

  int x_dims = x[0].size();
  int u_dims = u[0].size();
  VectorXd xt(x_dims);
  VectorXd ut(u_dims);

  //#pragma omp parallel for private(xt,ut)
  for(int t=start_T; t<end_T; t++) {
    xt = x[t];
    ut = (t<T) ? u[t] : VectorXd::Zero(u_dims);

    if(t<T) 
      c = [this, &ut](VectorXd xt){return model->cost(xt,ut);};
    else 
      c = [this](VectorXd xt){return model->final_cost(xt);};

    finite_diff_hessian(c, xt, c_xx[t]);
  }
}

void iLQR::calculate_cuu(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_uu, int start_T, int end_T) {
  int x_dims = x[0].size();
  int u_dims = u[0].size();
  VectorXd xt(x_dims);
  VectorXd ut(u_dims);

  //#pragma omp parallel for private(xt,ut)
  for(int t=start_T; t<end_T; t++) {
    xt = x[t];
    ut = (t<T) ? u[t] : VectorXd::Zero(u_dims);

    std::function<double(VectorXd)> c_t = [this, &xt](VectorXd ut){return model->cost(xt,ut);};
    finite_diff_hessian(c_t, ut, c_uu[t]);
  }
}

void iLQR::calculate_cxu(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xu, int start_T, int end_T) {
  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};
  std::function<double(VectorXd)> cf = std::bind(&Model::final_cost, model, _1);

  int x_dims = x[0].size();
  int u_dims = u[0].size();
  VectorXd xt(x_dims);
  VectorXd ut(u_dims);

  //#pragma omp parallel for private(ut)
  for(int t=start_T; t<end_T; t++) {
    xt = x[t];
    ut = (t<T) ? u[t] : VectorXd::Zero(u_dims);
    
    VectorXd px, mx, pu, mu;
    for (int i=0; i<x_dims; i++) {
      for (int j=0; j<u_dims; j++) {
        px = mx = x[t];
        pu = mu = ut;
        px(i) += eps2;
        mx(i) -= eps2;
        pu(j) += eps2;
        mu(j) -= eps2;
        if(t<T) 
          c_xu[t](i,j) = (c(px, pu) - c(mx, pu) - c(px, mu) + c(mx, mu)) / (4*sqr(eps2));
        else 
          c_xu[t](i,j) = (cf(px) - cf(mx) - cf(px) + cf(mx)) / (4*sqr(eps2)); // TODO this is wrong
      }
    }
  }
}
