#pragma once
#include "standardIncludes.h"
#define epsilon 1e-3

#define DISCRETE_TIME // this is equivalent to the discrete time formulation but f stays as: x_i+1=x_i+f(x_i,u_i)*dt rather than: x_i+1=f(x_i,u_i)
#define MITROVIC

#if defined(MITROVIC)
const double lambdaFactor = sqrt(10);
const int maxLambda = 1e7;
const double convergenceThreshold = 1e-15;
const int maxIterations = 200;
#else
const double lambdaFactor = 10;
const double convergenceThreshold = 1e-3;
const int maxLambda = 1000;
const int maxIterations = 100;
#endif
// Iterative Linear Quadratic Regulator
// variable names match the Linear Quadratic Regulator (wikipedia)
// which are also the Mitrovic 2010 (and Todorov 05) names unless stated otherwise
// Tassa variable names are the names given in Tassa/Todorov papers 2012,2014, and by StudyWolf 2016
// Equation numbers refer to the Tassa/Todorov 2012 paper
struct ILQR
{
  int n; // length of state
  int m; // length of signal
  int K; // trajectory length
  vector<VectorXd> xs;
  vector<VectorXd> us;
  vector<VectorXd> ls; // Tassa k
  vector<MatrixXd> Ls; // Tassa K
  VectorXd xd; // target state

  void generateFeedbackController(const VectorXd &x_0, const VectorXd &x_d, int trajectoryLength, int signalSize)
  {
    n = x_0.size();
    xd = x_d;
    K = trajectoryLength; // length of trajectory
    m = signalSize;
    xs.resize(K);
    us.resize(K);
    ls.resize(K);
    Ls.resize(K);
    for (int k = 0; k<K; k++)
    {
      xs[k] = VectorXd::Zero(n);
      us[k] = VectorXd::Zero(m);
      ls[k] = VectorXd::Zero(m);
      Ls[k] = MatrixXd::Zero(m,n);
    }
    xs[0] = x_0;
    double lambda = 1.0;
#if defined(MITROVIC)
    lambda = 100;
#endif
    bool simulateNewTrajectory = true;
    double cost = simulate(xs, us, ls, Ls, xs, us);
    vector<MatrixXd> A(K-1); // Tassa f_x
    vector<MatrixXd> B(K-1); // Tassa f_u

    vector<VectorXd> q(K); // Tassa l_x
    vector<VectorXd> r(K-1); // Tassa l_u
    vector<MatrixXd> Q(K); // Tassa l_xx
    vector<MatrixXd> R(K-1); // Tassa l_uu
    vector<MatrixXd> N(K-1); // Mitrovic P^T, Tassa L_ux^T
    double dt = timeDelta;
    for (int iterations = 0; iterations < maxIterations; iterations++)
    {
      // 1. Derivatives: Given a nominal (x,u,i) sequence, compute the derivatives of l and f in the RHS of eq (5)
      if (simulateNewTrajectory)
      {
        for (int k = 0; k<K-1; k++)
        {
          // x(t+1) = f(x(t), u(t)) = x(t) + dx(t) * dt
          // linearized dx(t) = np.dot(A(t), x(t)) + np.dot(B(t), u(t))
          getDynamicsDerivatives(xs[k], us[k], A[k], B[k]);
          getCostDerivatives(xs[k], us[k], q[k], r[k], Q[k], R[k], N[k]);
          // convert to discrete system
          A[k] = MatrixXd::Identity(n,n) + A[k]*dt;
          B[k] *= dt;
          q[k] *= dt; r[k] *= dt; Q[k] *= dt; R[k] *= dt; N[k] *= dt;
        }
        getFinalCostDerivatives(xs.back(), q.back(), Q.back());
        simulateNewTrajectory = false;
      }

      // 2. Backward pass: Iterate Eqs (5,10,11) for decreasing i=N-1,..1

      VectorXd p = q.back(); // Tassa: V_x.  NOTE: we're working backwards, so V_x = V_x[t+1] = V'_x
      MatrixXd P = Q.back(); // Tassa: V_xx
      for (int k = K-2; k>=0; k--)
      {
        MatrixXd H = R[k];
#if defined(DISCRETE_TIME)
        // Eq 5:
        H += B[k].transpose()*P*B[k];  // Tassa: Q_uu
#endif

        // Eq 10:
        // Calculate H^-1 with regularization term set by
        // Levenberg-Marquardt heuristic (at end of this loop)
        // [Note: this is different from the 2012 paper, more like 2014]
        SelfAdjointEigenSolver<MatrixXd> eigenSolver(H);
        ASSERT(eigenSolver.info() == Success);
        VectorXd eigenValues = eigenSolver.eigenvalues();
        for (int i = 0; i<eigenValues.size(); i++)
        {
          if (eigenValues[i] < 0.0)
            eigenValues[i] = 0.0;
          eigenValues[i] = 1.0 / (eigenValues[i] + lambda);
        }
        MatrixXd Hi = eigenSolver.eigenvectors() * eigenValues.asDiagonal() * eigenSolver.eigenvectors().transpose();

        // Eq 11:
        MatrixXd Hinv = Hi*H*Hi;
#if defined(MITROVIC) // also Todorov 2012 but not Tassa/Todorov 2014 (or studyWolf)
      //  Hinv = 2*Hi - Hinv; // corrected (assuming Hi is symmetric) // not apparenltly working
#endif
        // Eq 5:
        MatrixXd NBPA = N[k].transpose() + B[k].transpose()*P*A[k]; // Tassa: Q_ux, Mitrovic G
        Ls[k] = Hi*NBPA;
        P = Q[k] + A[k].transpose()*P*A[k] - NBPA.transpose()*Hinv*NBPA;

        VectorXd rBp = r[k] + B[k].transpose()*p; // Tassa: Q_u, Mitrovic g
        ls[k] = Hi*rBp;
        p = q[k] + A[k].transpose()*p - NBPA.transpose()*Hinv*rBp;
      }

      // 3. Forward pass: Set alpha=1. Iterate (12) and (8c) to compute a new nominal sequence.
      vector<VectorXd> xs2, us2;
      double costnew = simulate(xs, us, ls, Ls, xs2, us2);
      // If the integration diverged or condition (13) was not met, decrease alpha and restart the forward pass

      // Levenberg-Marquardt heuristic
      if (costnew < cost)
      {
        // decrease lambda (get closer to Newton's method)
        lambda /= lambdaFactor;

        xs = xs2; // update trajectory
        us = us2; // update control signal

        simulateNewTrajectory = true; // do another rollout

        // check to see if update is small enough to exit
        if (iterations > 0 && abs(costnew - cost) / cost < convergenceThreshold)
        {
          cost = costnew;
          cout << "Converged at iteration = " << iterations << "; Cost = " << cost << " logLambda = " << log(lambda) << endl;
          break;
        }
        cost = costnew;
      }
      else
      {
        // increase lambda (get closer to gradient descent)
        lambda *= lambdaFactor;
        if (lambda > maxLambda)
        {
          cout << "lambda > max_lambda at iteration = " << iterations << "; Cost = " << cost << "; logLambda = " << log(lambda) << endl;
          break;
        }
      }
    }
    cout << "Cost: " << cost << endl;
  }

  // OPTION A: overload these functions as needed with cost derivatives.
  // Note, in LQR syntax: f_x = A, f_u = B, l_xx = Q, l_uu = R, l_xu = N
  virtual void getDynamicsDerivatives(const VectorXd &x, const VectorXd &u, MatrixXd &f_x, MatrixXd &f_u)
  {
    f_x.resize(n,n);
    for (int i = 0; i<n; i++)
    {
      VectorXd eps = VectorXd::Zero(n);
      eps[i] = epsilon;
      f_x.col(i) = (f(x+eps, u) - f(x-eps, u)) / (2.0 * epsilon);
    }
    f_u.resize(n,m);
    for (int i = 0; i<m; i++)
    {
      VectorXd eps = VectorXd::Zero(m);
      eps[i] = epsilon;
      f_u.col(i) = (f(x, u+eps) - f(x, u-eps)) / (2.0 * epsilon);
    }
  }
  virtual void getCostDerivatives(const VectorXd &x, const VectorXd &u, VectorXd &l_x, VectorXd &l_u, MatrixXd &l_xx, MatrixXd &l_uu, MatrixXd &l_ux)
  {
    l_x.resize(n);
    l_u.resize(m);
    l_xx.resize(n,n);
    l_uu.resize(m,m);
    l_ux.resize(n,m);
    for (int i = 0; i<n; i++)
    {
      VectorXd eps = VectorXd::Zero(n);
      eps[i] = epsilon;
      l_x[i] = (l(x+eps, u) - l(x-eps, u)) / (2.0 * epsilon);
    }
    for (int i = 0; i<m; i++)
    {
      VectorXd eps = VectorXd::Zero(m);
      eps[i] = epsilon;
      l_u[i] = (l(x, u+eps) - l(x, u-eps)) / (2.0 * epsilon);
    }
    for (int i = 0; i<n; i++)
    {
      VectorXd eps = VectorXd::Zero(n);
      eps[i] = epsilon;
      for (int j = 0; j<n; j++) // ouch, this is the very slow part!
      {
        VectorXd eps2 = VectorXd::Zero(n);
        eps2[j] = epsilon;
        double d1 = (l(x+eps+eps2, u) - l(x-eps+eps2, u)) / (2.0 * epsilon);
        double d0 = (l(x+eps-eps2, u) - l(x-eps-eps2, u)) / (2.0 * epsilon);
        l_xx(j,i)= (d1-d0) / (2.0 * epsilon);
      }
    }
    for (int i = 0; i<m; i++)
    {
      VectorXd eps = VectorXd::Zero(m);
      eps[i] = epsilon;
      for (int j = 0; j<n; j++) // ouch, this is the very slow part!
      {
        VectorXd eps2 = VectorXd::Zero(n);
        eps2[j] = epsilon;
        double f1 = (l(x+eps2, u+eps) - l(x+eps2, u-eps)) / (2.0 * epsilon);
        double f0 = (l(x-eps2, u+eps) - l(x-eps2, u-eps)) / (2.0 * epsilon);
        l_ux(j,i)= (f1-f0) / (2.0 * epsilon);
      }
      for (int j = 0; j<m; j++) // ouch, this is the very slow part!
      {
        VectorXd eps2 = VectorXd::Zero(m);
        eps2[j] = epsilon;
        double e1 = (l(x, u+eps+eps2) - l(x, u-eps+eps2)) / (2.0 * epsilon);
        double e0 = (l(x, u+eps-eps2) - l(x, u-eps-eps2)) / (2.0 * epsilon);
        l_uu(j,i)= (e1-e0) / (2.0 * epsilon);
      }
    }
  }
  virtual void getFinalCostDerivatives(const VectorXd &x, VectorXd &l_x, MatrixXd &l_xx)
  {
    l_x.resize(n);
    l_xx.resize(n,n);
    for (int i = 0; i<n; i++)
    {
      VectorXd eps = VectorXd::Zero(n);
      eps[i] = epsilon;
      l_x[i] = (lf(x+eps) - lf(x-eps)) / (2.0 * epsilon);
    }
    for (int i = 0; i<n; i++)
    {
      VectorXd eps = VectorXd::Zero(n);
      eps[i] = epsilon;
      for (int j = 0; j<n; j++) // ouch, this is the very slow part!
      {
        VectorXd eps2 = VectorXd::Zero(n);
        eps2[j] = epsilon;
        double d1 = (lf(x+eps+eps2) - lf(x-eps+eps2)) / (2.0 * epsilon);
        double d0 = (lf(x+eps-eps2) - lf(x-eps-eps2)) / (2.0 * epsilon);
        l_xx(j,i)= (d1-d0) / (2.0 * epsilon);
      }
    }
  }

  // OPTION B: overload dynamics and cost functions as necessary.
  // Derivatives will be calculated by finite differences above
  // xdot = f(x,u)
  virtual VectorXd f(const VectorXd &x, const VectorXd &u) // non-linear state dynamics
  {
    // simple second order motion
    VectorXd xdot = VectorXd::Zero(n);
    xdot.head(n/2) = x.tail(n/2); // head is the position part of x, tail is the velocity part
    return xdot;
  }
  // total cost = lf(xf) + int l(x,u) dt
  virtual double l(const VectorXd &x, const VectorXd &u) // non-linear cost function
  {
    // just penalise the signal here
    const double forceCost = 1.0;
    return forceCost * u.transpose()*u;
  }
  virtual double lf(const VectorXd &x) // non-linear final cost
  {
    // penalise the end position and end velocity;
    const double positionCost = 1.0;
    const double velocityCost = 1.0;
    VectorXd pos = (xd-x).head(n/2);
    VectorXd vel = (xd-x).tail(n/2);
    return positionCost * pos.dot(pos) + velocityCost * vel.dot(vel);
  }

  // apply the feedback control to the whole trajectory
  double simulate(const vector<VectorXd> xs, const vector<VectorXd> us, const vector<VectorXd> &ls, const vector<MatrixXd> &Ls, vector<VectorXd> &xs2, vector<VectorXd> &us2)
  {
    xs2.resize(K);
    us2.resize(K);
    xs2[0] = xs[0];
    double cost = 0;

    // Run simulation with substeps
    for (int k = 0; k < K-1; k++)
    {
      // this is a relative LQR (so it operates on deltas from the nominal trajectory)
      VectorXd dx = xs2[k] - xs[k];
      VectorXd du = -ls[k] - Ls[k]*dx;

      us2[k] = us[k] + du;
      cost += l(xs2[k], us2[k])*timeDelta;
      xs2[k+1] = xs2[k] + f(xs2[k], us2[k])*timeDelta;
    }
    // Adjust for final cost, subsample trajectory
    cost += lf(xs2.back());

    return cost;
  }

  // this applies the feedback control one frame at a time, so external perturbations can be made to state
  void updateState(VectorXd &state, int k)
  {
    VectorXd dx = state - xs[k];
    VectorXd du = -ls[k] - Ls[k]*dx;
    state += f(state, us[k] + du)*timeDelta;
  }
};
