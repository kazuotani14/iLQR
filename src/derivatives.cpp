#include "ilqr.h"
#include "finite_diff.h"

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
    //TODO figure out how to avoid copying
    VectorXd xi = x[t];
    VectorXd ui = u[t];
    std::function<VectorXd(VectorXd)> dyn_x = [ui, this](VectorXd x_v){return model->dynamics(x_v,ui);};
    std::function<VectorXd(VectorXd)> dyn_u = [xi, this](VectorXd u_v){return model->dynamics(xi,u_v);};

    fx[t] = finite_diff_vec2vec(dyn_x, x[t], 4);
    fu[t] = finite_diff_vec2vec(dyn_u, u[t], 4);
  }
}

void iLQR::get_cost_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
															   VecOfVecXd &cx, VecOfVecXd &cu)
{
  for (int t=0; t<T; t++)
  {
    //TODO figure out how to avoid copying
    VectorXd xi = x[t];
    VectorXd ui = u[t];
    std::function<double(VectorXd)> cost_x = [ui, this](VectorXd x_v){return model->cost(x_v,ui);};
    std::function<double(VectorXd)> cost_u = [xi, this](VectorXd u_v){return model->cost(xi,u_v);};

    cx[t] = finite_diff_vec2scalar(cost_x, x[t]);
    cu[t] = finite_diff_vec2scalar(cost_u, u[t]);
  }
}

void iLQR::get_cost_2nd_derivatives(const VecOfVecXd &x, const VecOfVecXd &u,
                              VecOfMatXd &cxx, VecOfMatXd &cxu, VecOfMatXd &cuu)
{
  for (int t=0; t<T; t++)
  {
    VectorXd xi = x[t];
    VectorXd ui = u[t];



  }
}
//
//     //TODO remove repetition
//
//     VecXd pp, pm, mp, mm; //plus-plus, plus-minus, ....
//
//     //cxx
//     for (int i=0; i<n; i++){
//       for (int j=i; j<n; j++){
//         pp = pm = mp = mm = x[t];
//         pp(i) = pm(i) += eps;
//         mp(i) = mm(i) -= eps;
//         pp(j) = mp(j) += eps;
//         pm(j) = mm(j) -= eps;
//         cxx[t](i,j) = cxx[t](j,i) = (cost(pp, u[t]) - cost(mp, u[t]) - cost(pm, u[t]) + cost(mm, u[t])) / (4*sqr(eps));
//       }
//     }
//     //cxu
//     VecXd px, mx, pu, mu;
//     for (int i=0; i<n; i++){
//       for (int j=0; j<m; j++){
//         px = mx = x[t];
//         pu = mu = u[t];
//         px(i) += eps;
//         mx(i) -= eps;
//         pu(j) += eps;
//         mu(j) -= eps;
//         cxu[t](i,j) = (cost(px, pu) - cost(mx, pu) - cost(px, mu) + cost(mx, mu)) / (4*sqr(eps));
//       }
//     }
//     //cuu
//     for (int i=0; i<m; i++){
//       for (int j=i; j<m; j++){
//         pp = pm = mp = mm = u[t];
//         pp(i) = pm(i) += eps;
//         mp(i) = mm(i) -= eps;
//         pp(j) = mp(j) += eps;
//         pm(j) = mm(j) -= eps;
//         cuu[t](i,j) = cuu[t](j,i) = (cost(x[t], pp) - cost(x[t], mp) - cost(x[t], pm) + cost(x[t], mm)) / (4*sqr(eps));
//       }
//     }
//
//     // std::cout << t << '\n';
//     // std::cout << "cxx[t]: \n" << cxx[t] << '\n';
//     // std::cout << "cxu[t]: \n" << cxu[t] << '\n';
//     // std::cout << "cuu[t]: \n" << cuu[t] << '\n';
//     // getchar();
//   }
// } //get_cost_2nd_derivatives
