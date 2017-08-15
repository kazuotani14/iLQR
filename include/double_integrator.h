#ifndef _DOUBLE_INTEGRATOR_H_
#define _DOUBLE_INTEGRATOR_H_

#include "model.h"

class DoubleIntegrator : public Model
{

public:
  DoubleIntegrator(VectorXd xd): Hx(4,4), Hu(2,2), goal(xd)
  {
    x_dims = 4;
    u_dims = 2;

    Hx << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 0.2, 0,
         0, 0, 0, 0.2;
    Hu << 1, 0,
          0, 1;
    u_min = Vector2d(-0.5, -0.5);
    u_max = Vector2d(0.5, 0.5);
  }

  virtual VectorXd dynamics(const VectorXd& x, const VectorXd& u) override
  {
    VectorXd dx(4);
    dx(0) = x(2);
    dx(1) = x(3);
    dx(2) = u(0)/mass;
    dx(3) = u(1)/mass;

    return dx;
  }

  virtual double cost(const VectorXd& x, const VectorXd& u) override
  {
    double cost_x = (goal-x).transpose()*Hx*(goal-x);
    double cost_u = u.transpose()*Hu*u;
    return cost_x + cost_u;
  }

  virtual double final_cost(const VectorXd& x) override
  {
    double cost = (goal-x).transpose()*Hx*(goal-x);
    return cost;
  }

  virtual VectorXd integrate_dynamics(const VectorXd& x, const VectorXd& u, double dt) override
  {
    VectorXd x1 = x + dynamics(x,u)*dt;
    return x1;
  }


private:
  double mass = 1.0;
  MatrixXd Hx, Hu;
  VectorXd goal;
};

#endif
