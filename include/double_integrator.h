#ifndef _DOUBLE_INTEGRATOR_H_
#define _DOUBLE_INTEGRATOR_H_

#include "model.h"

class DoubleIntegrator : public Model
{

public:
    DoubleIntegrator(): Hx(4,4), Hu(2,2), goal(4)
    {
        x_dims = u_dims = 2;

        Hx << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 0.2, 0,
             0, 0, 0, 0.2;
        Hu << 1, 0,
              0, 1;
        goal << 1.0, 1.0, 0.0, 0.0;
        u_min = Vector2d(-0.5, -0.5);
        u_max = Vector2d(0.5, 0.5);
    }

    virtual VectorXd dynamics(const VectorXd& x, const VectorXd& u) override;
    virtual double cost(const VectorXd& x, const VectorXd& u) override;
    virtual double final_cost(const VectorXd& x) override;

    // this doesn't have to be virtual?
    virtual VectorXd integrate_dynamics(const VectorXd& x, const VectorXd& u, double dt) override;

private:
    double mass = 1.0;
    MatrixXd Hx, Hu;
    VectorXd goal;
};

#endif
