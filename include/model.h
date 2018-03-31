#ifndef _DYNAMICS_H_
#define _DYNAMICS_H_

#include "common.h"

class Model {
public:
  virtual VectorXd dynamics(const VectorXd& x, const VectorXd& u) = 0;
  virtual double cost(const VectorXd& x, const VectorXd& u) = 0;
  virtual double final_cost(const VectorXd& x) = 0;

  VectorXd integrate_dynamics(const VectorXd& x, const VectorXd& u, double dt) {
    VectorXd x1 = x + dynamics(x,u)*dt;
    return x1;
  }

  VectorXd u_min, u_max;

  int x_dims;
  int u_dims;
};

#endif
