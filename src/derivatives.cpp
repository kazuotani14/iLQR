#include "ilqr.h"
#include "finite_diff.h"
#include <thread>
#include <functional>

#define eps 1e-6 // For finite differencing
#define eps2 1e-3

using namespace std::placeholders;

// Note: this is not used in main ilqr anymore, mostly just here cuz tests use it
// Given a trajectory {x(t),u(t)} from forward pass, compute deriviatives along it
void iLQR::compute_derivatives(const VecOfVecXd& x, const VecOfVecXd& u)
{
  get_dynamics_derivatives(x, u, fx, fu);
  get_cost_derivatives(x, u, cx, cu);
  // get_cost_2nd_derivatives(x, u, cxx, cxu, cuu);
  get_cost_2nd_derivatives_mt(x, u, cxx, cxu, cuu, 10);
}

// Updates fx, fu
void iLQR::get_dynamics_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& f_x, VecOfMatXd& f_u)
{
  for (int t=0; t<T; t++)
  {
    //TODO figure out how to avoid copying
    VectorXd xi = x[t];
    VectorXd ui = u[t];

    std::function<VectorXd(VectorXd)> dyn_x = std::bind(&Model::integrate_dynamics, model, _1, u[t], dt);
    std::function<VectorXd(VectorXd)> dyn_u = std::bind(&Model::integrate_dynamics, model, x[t], _1, dt);

    f_x[t] = finite_diff_vec2vec(dyn_x, x[t], 4);
    f_u[t] = finite_diff_vec2vec(dyn_u, u[t], 4);
  }
}

// Updates cx, cu
void iLQR::get_cost_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfVecXd& c_x, VecOfVecXd& c_u)
{
  for (int t=0; t<T+1; t++)
  {
    Vector2d ut;
    if(t<T)
    {
      ut = u[t];
    }
    else
    {
      ut = Vector2d::Zero();
    }
    std::function<double(VectorXd)> cost_x = std::bind(&Model::cost, model, _1, ut);
    std::function<double(VectorXd)> cost_u = std::bind(&Model::cost, model, x[t], _1);
    std::function<double(VectorXd)> cost_f = std::bind(&Model::final_cost, model, _1);

    if(t<T)
    {
      c_x[t] = finite_diff_vec2scalar(cost_x, x[t]);
      c_u[t] = finite_diff_vec2scalar(cost_u, ut);
    }
    else
    {
      c_x[t] = finite_diff_vec2scalar(cost_f, x[t]);
      c_u[t].resize(2);
      c_u[t].setZero();
    }
  }
}

// Update cxx, cxu, cuu
void iLQR::get_cost_2nd_derivatives(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xx, VecOfMatXd& c_xu, VecOfMatXd& c_uu)
{
  VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  int n = x[0].size();
  int m = u[0].size();
  for (int t=0; t<T+1; t++)
  {
    c_xx[t].resize(n,n);
    c_xu[t].resize(n,m);
    c_uu[t].resize(m,m);
  }

  calculate_cxx(x, u, c_xx, 0, T+1);
  calculate_cxu(x, u, c_xu, 0, T+1);
  calculate_cuu(x, u, c_uu, 0, T+1);
}

void iLQR::get_cost_2nd_derivatives_mt(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xx, VecOfMatXd& c_xu, VecOfMatXd& c_uu, int n_threads_per)
{
  assert((T+1)%n_threads_per == 0);

  //Setup
  int n = x[0].size();
  int m = u[0].size();
  for (int t=0; t<T+1; t++)
  {
    c_xx[t].resize(n,n);
    c_xu[t].resize(n,m);
    c_uu[t].resize(m,m);
  }

  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  //Spawn worker threads
  std::vector<std::thread> threads;
  for(int i=0; i<n_threads_per; i++)
  {
    int start_T = (T+1)/n_threads_per * i;
    int end_T = (T+1)/n_threads_per * (i+1);

    threads.emplace_back(std::bind(&iLQR::calculate_cxx, this, x, u, c_xx, start_T, end_T));
    threads.emplace_back(std::bind(&iLQR::calculate_cxu, this, x, u, c_xu, start_T, end_T));
    threads.emplace_back(std::bind(&iLQR::calculate_cuu, this, x, u, c_uu, start_T, end_T));
  }

  for(std::thread& t : threads) t.join();
}

void iLQR::calculate_cxx(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xx, int start_T, int end_T)
{
  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

  int n = x[0].size();
  for(int t=start_T; t<end_T; t++)
  {
    Vector2d ut;
    if(t==T) ut << 0., 0.;
    else ut = u[t];

    for (int i=0; i<n; i++){
      for (int j=i; j<n; j++){
        pp = pm = mp = mm = x[t];
        pp(i) += eps2;
        pp(j) += eps2;
        pm(i) += eps2;
        pm(j) -= eps2;
        mp(i) -= eps2;
        mp(j) += eps2;
        mm(i) -= eps2;
        mm(j) -= eps2;
        c_xx[t](i,j) = cxx[t](j,i) = (c(pp, ut) - c(mp, ut) - c(pm, ut) + c(mm, ut)) / (4*sqr(eps2));
      }
    }
  }
}

void iLQR::calculate_cxu(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_xu, int start_T, int end_T)
{
  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  VectorXd px, mx, pu, mu;

  int n = x[0].size();
  int m = u[0].size();
  for(int t=start_T; t<end_T; t++)
  {
    Vector2d ut;
    if(t==T) ut << 0., 0.;
    else ut = u[t];

    for (int i=0; i<n; i++){
      for (int j=0; j<m; j++){
        px = mx = x[t];
        pu = mu = ut;
        px(i) += eps2;
        mx(i) -= eps2;
        pu(j) += eps2;
        mu(j) -= eps2;
        c_xu[t](i,j) = (c(px, pu) - c(mx, pu) - c(px, mu) + c(mx, mu)) / (4*sqr(eps2));
      }
    }
  }
}

void iLQR::calculate_cuu(const VecOfVecXd& x, const VecOfVecXd& u, VecOfMatXd& c_uu, int start_T, int end_T)
{
  std::function<double(VectorXd, VectorXd)> c = [this](VectorXd x_v, VectorXd u_v){return model->cost(x_v,u_v);};

  VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....

  int m = u[0].size();
  for(int t=start_T; t<end_T; t++)
  {
    Vector2d ut;
    if(t==T) ut << 0., 0.;
    else ut = u[t];

    for (int i=0; i<m; i++){
      for (int j=i; j<m; j++){
        pp = pm = mp = mm = ut;
        pp(i) += eps2;
        pp(j) += eps2;
        pm(i) += eps2;
        pm(j) -= eps2;
        mp(i) -= eps2;
        mp(j) += eps2;
        mm(i) -= eps2;
        mm(j) -= eps2;
        c_uu[t](i,j) = cuu[t](j,i) = (c(x[t], pp) - c(x[t], mp) - c(x[t], pm) + c(x[t], mm)) / (4*sqr(eps2));
      }
    }
  }
}
