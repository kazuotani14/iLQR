#include "ilqr.h"
#include "double_integrator.h"

#define SHOWPROGRESS
// #define VERBOSE

double iLQR::init_traj(const VectorXd &x_0, const VecOfVecXd &u_0)
{
  T = u_0.size();

  //initialize xs and us
  xs.resize(T+1);
  us.resize(T);
  x0 = x_0;

  // call forward_pass to get xs, us, cost
  double cost_i = forward_pass(x0, u_0);
  std::cout << "Initial cost: " << cost_i << "\n";
  cost_s = cost_i;

  //allocate space for later
  du = MatrixXd(2,T);
  fx.resize(T+1);
  fu.resize(T+1);
  cx.resize(T+1);
  cu.resize(T+1);
  cxu.resize(T+1);
  cxx.resize(T+1);
  cuu.resize(T+1);

   dV = Vector2d(2,1);
  Vx.resize(T+1);
  Vxx.resize(T+1);
  k.resize(T);
  K.resize(T);

  int n = model->x_dims;
  int m = model->u_dims;

  std::fill(fx.begin(), fx.end(), MatrixXd::Zero(n,n));
  std::fill(fu.begin(), fu.end(), MatrixXd::Zero(n,m));
  std::fill(cx.begin(), cx.end(), VectorXd::Zero(n));
  std::fill(cu.begin(), cu.end(), VectorXd::Zero(m));
  std::fill(cxx.begin(), cxx.end(), MatrixXd::Zero(n,n));
  std::fill(cxu.begin(), cxu.end(), MatrixXd::Zero(n,m));
  std::fill(cuu.begin(), cuu.end(), MatrixXd::Zero(m,m));
  std::fill(Vx.begin(), Vx.end(), VectorXd::Zero(n));
  std::fill(Vxx.begin(), Vxx.end(), MatrixXd::Zero(n,n));
  std::fill(k.begin(), k.end(), VectorXd::Zero(m));
  std::fill(K.begin(), K.end(), MatrixXd::Zero(m,n));

  return cost_s;
}

// Initialize trajectory with control sequence
void iLQR::generate_trajectory(const VectorXd &x_0, const VecOfVecXd &u0)
{
  init_traj(x_0, u0);
  generate_trajectory();
}

// Warm-start
void iLQR::generate_trajectory(const VectorXd &x_0)
{
  assert(us.size() > 0);

  x0 = x_0;

  double cost_i = forward_pass(x_0, us);
  std::cout << "Initial cost: " << cost_i << "\n";
  cost_s = cost_i;

  generate_trajectory();
}

// This assumes that x0, xs, us are initialized
void iLQR::generate_trajectory()
{
  std::cout << "assertions" << std::endl;
  assert(x0.size()>0);
  assert(!xs.empty());
  assert(!us.empty());
  std::cout << "..." << std::endl;

  VecOfVecXd x_old, u_old;

  // constants, timers, counters
  bool flgChange = true;
  bool stop = false;
  double dcost = 0;
  double z = 0;
  double expected = 0;
  int diverge = 0;
  double new_cost, gnorm;

  #ifdef VERBOSE
    std::cout << "\n=========== begin iLQG ===========\n";
  #endif

  int iter;
  for (iter=0; iter<maxIter; iter++)
  {
    x_old = xs; u_old = us;

    if (stop)
      break;

    #ifdef VERBOSE
      std::cout << "Iteration " << iter << ".\n";
    #endif

    //--------------------------------------------------------------------------
    //STEP 1: Differentiate dynamics and cost along new trajectory

    // auto start = std::chrono::system_clock::now();

    if (flgChange)
    {
      get_dynamics_derivatives(xs, us, fx, fu);
      get_cost_derivatives(xs, us, cx, cu);
      // get_cost_2nd_derivatives(x, u, cxx, cxu, cuu);
      get_cost_2nd_derivatives_mt(xs, us, cxx, cxu, cuu, 10);
      flgChange = 0;
    }
    #ifdef VERBOSE
      cout << "Finished step 1 : compute derivatives." << endl;;
    #endif

    // auto now = std::chrono::system_clock::now();
    // long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    // cout << "compute_derivatives took: " << elapsed/1000. << " seconds." << endl;

    //--------------------------------------------------------------------------
    // STEP 2: Backward pass, compute optimal control law and cost-to-go

    bool backPassDone = false;
    while (!backPassDone)
    {
       // update Vx, Vxx, l, L, dV with back_pass
      diverge = 0;
      diverge = backward_pass();

      if (diverge != 0)
      {
        #ifdef VERBOSE
          std::cout << "Backpass failed at timestep " << diverge << ".\n";
        #endif

        dlambda   = std::max(dlambda * lambdaFactor, lambdaFactor);
        lambda    = std::max(lambda * dlambda, lambdaMin);
        if (lambda > lambdaMax)
            break;
        continue;
      }
      backPassDone = true;
    }

    // check for termination due to small gradient
    gnorm = get_gradient_norm(k, us);
    if (gnorm < tolGrad && lambda < 1e-5)
    {
      #ifdef SHOWPROGRESS
        cout << "\nSUCCESS: gradient norm < tolGrad\n" << endl;
      #endif
      break;
    }

    #ifdef VERBOSE
      std::cout << "Finished step 2 : backward pass. \n";
    #endif

    //--------------------------------------------------------------------------
    // STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost

    bool fwdPassDone = 0;
    VecOfVecXd xnew(T+1);
    VecOfVecXd unew(T);
    double alpha;

    if (backPassDone) //  serial backtracking line-search
    {
      for (int i=0; i<Alpha.size(); i++)
      {
        alpha = Alpha(i);
        VecOfVecXd u_plus_feedforward = add_bias_to_u(us, k, alpha);

        new_cost = forward_pass(x0, u_plus_feedforward);
        dcost    = cost_s - new_cost;
        expected = -alpha * (dV(0) + alpha*dV(1));

        if (expected>0)
        {
          z = dcost/expected;
        }
        else
        {
          z = sgn(dcost);
          #ifdef VERBOSE
            cout << "Warning: non-positive expected reduction: should not occur" << endl;
          #endif
        }

        if(z > zMin)
        {
          fwdPassDone = 1;
          break;
        }
      }

      if (!fwdPassDone)
      {
        #ifdef VERBOSE
          cout << "Forward pass failed" << endl;
        #endif
        alpha = 0.0; // signals failure of forward pass
      }
    }

    #ifdef VERBOSE
      cout << "Finished step 3 : forward pass." <<endl;
    #endif

    //--------------------------------------------------------------------------
    // STEP 4: accept step (or not), print status
    #ifdef SHOWPROGRESS
      if (iter==0)
       std::cout << "iteration\tcost\t\treduction\texpect\t\tgrad\t\tlog10(lambda)\n";
    #endif

     if (fwdPassDone)
     {
      #ifdef SHOWPROGRESS
        printf("%-12d\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                  iter, new_cost, dcost, expected, gnorm, log10(lambda));
      #endif

       // decrease lambda
       dlambda   = std::min(dlambda / lambdaFactor, 1/lambdaFactor);
       lambda    = lambda * dlambda * (lambda > lambdaMin);

       // accept changes
      // cout << "accepting new cost: " << new_cost << endl;
       cost_s          = new_cost;
       flgChange       = true;

       //terminate?
       if (dcost < tolFun)
      {
        #ifdef SHOWPROGRESS
           std::cout << "\nSUCCESS: cost change < tolFun\n";
        #endif
         break;
       }
     }
     else // no cost improvement
     {
      xs = x_old;
      us = u_old;

       // increase lambda
       dlambda  = std::max(dlambda * lambdaFactor, lambdaFactor);
       lambda   = std::max(lambda * dlambda, lambdaMin);

      #ifdef SHOWPROGRESS
        printf("%-12d\t%-12s\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                  iter,"NO STEP", dcost, expected, gnorm, log10(lambda));
      #endif

       // terminate?
       if (lambda > lambdaMax){
        #ifdef SHOWPROGRESS
           std::cout << "\nEXIT: lambda > lambdaMax\n";
        #endif
         break;
       }
    }

    #ifdef SHOWPROGRESS
      if (iter==maxIter) std::cout << "\nEXIT: Maximum iterations reached.\n";
    #endif

  } // end top-level for-loop

  output_to_csv("ilqr_result.csv");

} //generate_trajectory


// Saves new state and control sequence in xs, us
//TODO make inputs/outputs explicit here?
double iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u)
{
  double total_cost = 0;

  VectorXd x_curr = x0;
  VectorXd x1;
  VectorXd u_curr;

  VecOfVecXd x_new(T+1);
  x_new[0] = x0;

  for(int t=0; t<T; t++)
  {
    u_curr = u[t];

    if (K.size()>0)
    {
      VectorXd dx = x_new[t] - xs[t];
      u_curr += K[t]*dx; //apply LQR control gains
    }

    us[t] = clamp_to_limits(u_curr, model->u_min, model->u_max);
    x1 = model->integrate_dynamics(x_curr, u_curr, dt);
    total_cost += model->cost(x_curr,u_curr);

    x_new[t+1] = x1;
    x_curr = x1;
  }

  xs = x_new;
  total_cost += model->final_cost(xs[T]);

  return total_cost;
}


/*
   INPUTS
      cx: 2x(T+1)          cu: 2x(T+1)
      cuu: nxnx(T+1)        cxx: nxnx(T+1)  cuu: 2x2x(T+1)
      fx: nxnx(T+1)        fu: nx2x(T+1)    fxx: none
      fxu: None            fuu: none        u: 2xT
    OUTPUTS
      Vx: nx(T+1)      Vxx: nxnx(T+1)      k:mxT
      K: mxnxT         dV: 2x1
      diverge - returns 0 if it doesn't diverge, timestep where it diverges otherwise
*/
//TODO make inputs/outputs explicit here?
int iLQR::backward_pass()
{
  int n = model->x_dims;
  int m = model->u_dims;

  //cost-to-go at end
  Vx[T] = cx[T];
  Vxx[T] = cxx[T];

  VectorXd Qx(n), Qu(m);
  MatrixXd Qxx(n,n), Qux(m,n), Quu(m,m);
  VectorXd k_i(m);
  MatrixXd K_i(m,n);
  dV.setZero();

  for (int i=(T-1); i>=0; i--) // back up from end of trajectory
  {
    Qx  = cx[i] + (fx[i].transpose() * Vx[i+1]);
    Qu  = cu[i] + (fu[i].transpose() * Vx[i+1]);
    Qxx = cxx[i] + (fx[i].transpose() * Vxx[i+1] * fx[i]);
    Qux = cxu[i].transpose() + (fu[i].transpose() * Vxx[i+1] * fx[i]);
      Quu = cuu[i] + (fu[i].transpose() * Vxx[i+1] * fu[i]);

      MatrixXd Vxx_reg = Vxx[i+1];
    MatrixXd Qux_reg = cxu[i].transpose() + (fu[i].transpose() * Vxx_reg * fx[i]);
    MatrixXd QuuF = cuu[i] + (fu[i].transpose() * Vxx_reg * fu[i]) + (lambda*MatrixXd::Identity(2,2));

    VectorXd lower = model->u_min - us[i];
    VectorXd upper = model->u_max - us[i];

    boxQPResult res = boxQP(QuuF, Qu, k[std::min(i+1,T-1)], lower, upper);

    int result = res.result;
    k_i = res.x_opt;
    MatrixXd R = res.H_free;
    VectorXd v_free = res.v_free;

    if(result < 1) return i;

    K_i.setZero();
    if (v_free.any())
    {
      MatrixXd Lfree;
      Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, v_free));

      int row_i = 0;
      for(int k=0; k<v_free.size(); k++)
      {
        if(v_free(k))
          K_i.row(k) = Lfree.row(row_i++);
      }
    }

    // update cost-to-go approximation
    dV(0) += k_i.transpose()*Qu;
    dV(1) += 0.5*k_i.transpose()*Quu*k_i;

    Vx[i]  = Qx  + K_i.transpose()*Quu*k_i + K_i.transpose()*Qu + Qux.transpose()*k_i;
    Vxx[i] = Qxx + K_i.transpose()*Quu*K_i + K_i.transpose()*Qux + Qux.transpose()*K_i;
    Vxx[i] = 0.5 * (Vxx[i] + Vxx[i].transpose());

      // save controls/gains
      k[i]     = k_i;
      K[i]     = K_i;
  }

  return 0;
}


VecOfVecXd iLQR::add_bias_to_u(const VecOfVecXd &u, const VecOfVecXd &l, const double alpha)
{
  VecOfVecXd new_u = u;
  for(unsigned int i=0; i<u.size(); i++){
    new_u[i] += l[i]*alpha;
  }
  return new_u;
}


// Replaces this line from matlab: g_norm = mean(max(abs(l) ./ (abs(u)+1), [], 1));
double iLQR::get_gradient_norm(const VecOfVecXd& l, const VecOfVecXd& u)
{
  std::vector<double> vals(l.size());
  for (unsigned int i=0; i<l.size(); i++)
  {
    VectorXd v = l[i].cwiseAbs().array() / (u[i].cwiseAbs().array() + 1);
    double max_val = v.maxCoeff();
    vals[i] = max_val;
  }
  double average = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
  return average;
}


void iLQR::output_to_csv(const std::string filename)
{
  FILE *XU = fopen(filename.c_str(), "w");

  for(unsigned int i=1; i<=xs[0].size(); i++)
    fprintf(XU, "x%d, ", i);
  for(unsigned int j=0; j<us[0].size(); j++)
    fprintf(XU, "u%d, ", j);
  fprintf(XU, "u%d\n", int(us[0].size()));

  for(int t=0; t<T; t++)
  {
    for(unsigned int i=0; i<xs[t].size(); i++)
      fprintf(XU, "%f, ", xs[t](i));
    for(unsigned int j=0; j<us[t].size()-1; j++)
      fprintf(XU, "%f, ", us[t](j));
    fprintf(XU, "%f\n", us[t](us[t].size()-1));
  }

  for(unsigned int i=0; i<xs[T].size(); i++)
    fprintf(XU, "%f, ", xs[T](i));

  fclose(XU);
  cout << "Saved iLQR result to " << filename << endl;
}
