#include "double_integrator.h"

VectorXd DoubleIntegrator::dynamics(const VectorXd& x, const VectorXd& u)
{
    VectorXd dx(4);
    dx(0) = x(2);
    dx(1) = x(3);
    dx(2) = u(0)/mass;
    dx(3) = u(1)/mass;

    return dx;
}

double DoubleIntegrator::cost(const VectorXd& x, const VectorXd& u)
{
    double cost_x = (goal-x).transpose()*Hx*(goal-x);
    double cost_u = u.transpose()*Hu*u;
    return cost_x + cost_u;
}

double DoubleIntegrator::final_cost(const VectorXd& x)
{
    double cost = x.transpose()*Hx*x;
    return cost;
}

VectorXd DoubleIntegrator::integrate_dynamics(const VectorXd& x, const VectorXd& u, double dt)
{
    VectorXd x1 = x + dynamics(x,u)*dt;
    return x1;
}
