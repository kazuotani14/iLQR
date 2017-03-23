/*
Basically just need to provide class with f(dynamics), l(cost), lf(final cost)
*/

#pragma once
#include "standardIncludes.h"
#include "iLQR.h"

class Acrobot : public iLQR
{
  double I1, I2;
  double l1,l2,m1,m2,g;
  double lc1, lc2;
public:
  void init()
  {
    g = 9.8;
    I1=I2=l1=l2=m1=m2=1;
    // derived quantities
    lc1 = 0.5*l1;
    lc2 = 0.5*l2;
  }
  // lagrangian dynamics:
  Matrix2d H(Vector2d &q)
  {
    double c2 = cos(q[1]);
    Matrix2d H;
    H(0,0) = I1+I2+m2*l1*l1+2*m2*l1*lc2*c2;
    H(0,1) = I2+m2*l1*lc2*c2;
    H(1,0) = I2+m2*l1*lc2*c2;
    H(1,1) = I2;
    return H;
  }
  Matrix2d C(Vector2d &q, Vector2d &qdot)
  {
    Matrix2d C;
    double s2 = sin(q[1]);
    C(0,0) = -2*m2*l1*lc2*s2*qdot[1];
    C(0,1) = -m2*l2*lc2*s2*qdot[1];
    C(1,0) = m2*l1*lc2*s2*qdot[0];
    C(1,1) = 0;
    return C;
  }
  Vector2d G(Vector2d &q)
  {
    Vector2d G;
    double s1 = sin(q[0]);
    double s1p2 = sin(q[0] + q[1]);
    G[0] = m1*g*lc1*s1 + m2*g*(l1*s1+lc2*s1p2);
    G[1] = m2*g*lc2*s1p2;
    return G;
  }

  VectorXd f(const VectorXd &x, const VectorXd &u) // non-linear state dynamics
  {
    // simple second order motion
    Vector2d q = x.head(2);
    Vector2d qdot = x.tail(2);
    VectorXd xdot(4);
    xdot.head(2) = qdot; // velocity = velocity
//    MatrixXd M = H(q).inverse(); // acceleration = m^-1 * (centrifugal/coriolis forces)
//    VectorXd u2 =  - C(q,qdot)*qdot; // acceleration = m^-1 * (centrifugal/coriolis forces)
//    VectorXd g2 = G(q); // acceleration = m^-1 * (centrifugal/coriolis forces)
    xdot.tail(2) = H(q).inverse()*(Vector2d(0,u[0]) - C(q,qdot)*qdot-G(q)); // acceleration = m^-1 * (centrifugal/coriolis forces)
    return xdot;
  }
  // total cost = lf(xf) + int l(x,u) dt
  virtual double l(const VectorXd &x, const VectorXd &u) // non-linear cost function
  {
    // not this is linear so we could implement it more efficiently by deriving the getCostDerivatives function directly
    // but this way it is easy to make it non-linear if needed
    Vector2d q = (xd-x).head(2);
    Vector2d qdot = (xd-x).tail(2);
    double Ks = 0.0; // penalise error
    double Kd = 0.0; // penalise velocity
    double Kr = 1.0; // penalise torque
    return Ks*Ks*q.dot(q) + Kd*Kd*qdot.dot(qdot) + Kr*Kr*u.dot(u);
  }

  virtual double lf(const VectorXd &x)
  {
    Vector2d q = (xd-x).head(2);
    Vector2d qdot = (xd-x).tail(2);
    double Ks = 20.0; // penalise error
    double Kd = 20.0; // penalise velocity
    return Ks*Ks*q.dot(q) + Kd*Kd*qdot.dot(qdot);
  }

  vector<Vector2d> getNodePositions(const Vector2d &q)
  {
    vector<Vector2d> pos(3);
    pos[0] << 0,0;
    pos[1] << l1*sin(q[0]), -l1*cos(q[0]);
    pos[2] << l2*sin(q[0]+q[1]), -l2*cos(q[0]+q[1]);
    pos[2] += pos[1];
    return pos;
  }
};
