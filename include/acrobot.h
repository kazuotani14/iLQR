#ifndef _ACROBOT_H_
#define _ACROBOT_H_

#include "model.h"

using Eigen::Matrix2d;
typedef Eigen::Matrix<double, 1, 1> Vector1d;

/*
  state = [theta1, theta2, dtheta1, dtheta2]
  theta1 defined wrt axis pointing down (balancing up is theta1=pi)
  theta2 defined wrt first link

  Borrowed from: https://github.com/TGlad/ILQR/blob/master/acrobot.h
*/
class Acrobot : public Model {

public:
  Acrobot(): g(9.81), Hx(4,4), Hu(2,2) {
    goal = VectorXd(4);
    goal << 3.1415, 0.0, 0.0, 0.0;

    I1=I2=l1=l2=m1=m2 = 1;
    lc1 = 0.5*l1;
    lc2 = 0.5*l2;

    x_dims = 4;
    u_dims = 1;

    Hx << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 0.2, 0,
          0, 0, 0, 0.2;
    Hu << 1, 0,
          0, 1;
    u_min = Vector1d(-100);
    u_max = Vector1d(100);
  }

  // lagrangian dynamics:
  Matrix2d H(Vector2d &q) {
    double c2 = cos(q[1]);
    Matrix2d H;
    H(0,0) = I1+I2+m2*l1*l1+2*m2*l1*lc2*c2;
    H(0,1) = I2+m2*l1*lc2*c2;
    H(1,0) = I2+m2*l1*lc2*c2;
    H(1,1) = I2;
    return H;
  }

  Matrix2d C(Vector2d &q, Vector2d &qdot) {
    Matrix2d C;
    double s2 = sin(q[1]);
    C(0,0) = -2*m2*l1*lc2*s2*qdot[1];
    C(0,1) = -m2*l2*lc2*s2*qdot[1];
    C(1,0) = m2*l1*lc2*s2*qdot[0];
    C(1,1) = 0;
    return C;
  }

  Vector2d G(Vector2d &q) {
    Vector2d G;
    double s1 = sin(q[0]);
    double s1p2 = sin(q[0] + q[1]);
    G[0] = m1*g*lc1*s1 + m2*g*(l1*s1+lc2*s1p2);
    G[1] = m2*g*lc2*s1p2;
    return G;
  }

  virtual VectorXd dynamics(const VectorXd& x, const VectorXd& u) override {
    // simple second order motion
    VectorXd dx(4);
    Vector2d q = x.head(2);
    Vector2d qdot = x.tail(2);

    dx.head(2) = qdot; // velocity = velocity
    dx.tail(2) = H(q).inverse()*(Vector2d(0,u[0]) - C(q,qdot)*qdot-G(q)); // acceleration = m^-1 * (centrifugal/coriolis forces)
    return dx;
  }

  virtual double cost(const VectorXd& x, const VectorXd& u) override {
    // quadratic costs
    Vector2d q = (goal-x).head(2);
    Vector2d qdot = (goal-x).tail(2);
    double Ks = 0.0; // penalize error
    double Kd = 0.0; // penalize velocity
    double Kr = 0.1; // penalize torque
    double cost_val = Ks*Ks*q.dot(q) + Kd*Kd*qdot.dot(qdot) + Kr*Kr*u.dot(u); 
    return cost_val;
  }

  virtual double final_cost(const VectorXd& x) override {
    Vector2d q = (goal-x).head(2);
    Vector2d qdot = (goal-x).tail(2);
    double Ks = 20.0; // penalize error
    double Kd = 20.0; // penalize velocity
    return Ks*Ks*q.dot(q) + Kd*Kd*qdot.dot(qdot); 
  }

  virtual VectorXd integrate_dynamics(const VectorXd& x, const VectorXd& u, double dt) override {
    VectorXd x1 = x + dynamics(x,u)*dt;
    return x1;
  }

  // vector<Vector2d> getNodePositions(const Vector2d &q)
  // {
  //   vector<Vector2d> pos(3);
  //   pos[0] << 0,0;
  //   pos[1] << l1*sin(q[0]), -l1*cos(q[0]);
  //   pos[2] << l2*sin(q[0]+q[1]), -l2*cos(q[0]+q[1]);
  //   pos[2] += pos[1];
  //   return pos;
  // }

private:
  double I1, I2;
  double l1, l2, m1, m2, g;
  double lc1, lc2;

  MatrixXd Hx, Hu; //LQR cost matrices
  VectorXd goal; 
};

#endif
