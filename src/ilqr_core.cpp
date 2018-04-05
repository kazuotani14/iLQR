#define SHOWPROGRESS
#define TIMESTUFF

#include "ilqr.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

double iLQR::init_traj(const VectorXd& x_0, const VecOfVecXd& u_0) {
  T = u_0.size();

  //initialize xs and us
  xs.resize(T+1);
  x0 = xs[0] = x_0;
  us = u_0;

  // call forward_pass to get xs, us, cost
  double cost_i = forward_pass(x0, us);

  //allocate everything
  du = MatrixXd(model->u_dims,T);
  fx.resize(T+1);
  fu.resize(T+1);
  cx.resize(T+1);
  cu.resize(T+1);
  cxu.resize(T+1);
  cxx.resize(T+1);
  cuu.resize(T+1);

  dV = Vector2d(model->u_dims,1);
  Vx.resize(T+1);
  Vxx.resize(T+1);
  k.resize(T);
  K.resize(T);

  std::fill(fx.begin(), fx.end(), MatrixXd::Zero(model->x_dims,model->x_dims));
  std::fill(fu.begin(), fu.end(), MatrixXd::Zero(model->x_dims,model->u_dims));
  std::fill(cx.begin(), cx.end(), VectorXd::Zero(model->x_dims));
  std::fill(cu.begin(), cu.end(), VectorXd::Zero(model->u_dims));
  std::fill(cxx.begin(), cxx.end(), MatrixXd::Zero(model->x_dims,model->x_dims));
  std::fill(cxu.begin(), cxu.end(), MatrixXd::Zero(model->x_dims,model->u_dims));
  std::fill(cuu.begin(), cuu.end(), MatrixXd::Zero(model->u_dims,model->u_dims));
  std::fill(Vx.begin(), Vx.end(), VectorXd::Zero(model->x_dims));
  std::fill(Vxx.begin(), Vxx.end(), MatrixXd::Zero(model->x_dims,model->x_dims));
  std::fill(k.begin(), k.end(), VectorXd::Zero(model->u_dims));
  std::fill(K.begin(), K.end(), MatrixXd::Zero(model->u_dims,model->x_dims));

#ifdef SHOWPROGRESS
  cout << "Initial cost: " << cost_i << endl;
#endif
  cost_s = cost_i;

  return cost_s;
}

// Initialize trajectory with control sequence
void iLQR::generate_trajectory(const VectorXd &x_0, const VecOfVecXd &u0) {
  init_traj(x_0, u0);
  generate_trajectory();
}

// Warm-start
void iLQR::generate_trajectory(const VectorXd &x_0) {
  assert(us.size() > 0);
  x0 = x_0;

  double cost_i = forward_pass(x_0, us);
#ifdef SHOWPROGRESS
  std::cout << "Initial cost: " << cost_i << std::endl;;
#endif
  cost_s = cost_i;

  generate_trajectory();
}

// This assumes that x0, xs, us are initialized
void iLQR::generate_trajectory() {  
  assert(x0.size()>0);
  assert(!xs.empty());
  assert(!us.empty());

  VecOfVecXd x_old, u_old;

  // constants, timers, counters
  bool flgChange = true;
  bool stop = false;
  double dcost = 0;
  double z = 0;
  double expected = 0;
  int diverge = 0;
  double new_cost, gnorm;

  #ifdef TIMESTUFF
    double t_compute_deriv = 0;
    double t_backward = 0;
    double t_forward = 0;
    auto all_start = std::chrono::system_clock::now();
  #endif

  int iter;
  for (iter=0; iter<maxIter; iter++) {
    x_old = xs; u_old = us;

    if (stop) break;

    //--------------------------------------------------------------------------
    //STEP 1: Differentiate dynamics and cost along new trajectory

    #ifdef TIMESTUFF
    auto start = std::chrono::system_clock::now();
    #endif

    if (flgChange) {
      get_dynamics_derivatives(xs, us, fx, fu);
      get_cost_derivatives(xs, us, cx, cu);
      get_cost_2nd_derivatives(xs, us, cxx, cxu, cuu);
      // get_cost_2nd_derivatives_mt(xs, us, cxx, cxu, cuu, 10);
      flgChange = 0;
    }

    #ifdef TIMESTUFF
    auto now = std::chrono::system_clock::now();
    long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    // cout << "compute_derivatives took: " << elapsed/1000. << " seconds." << endl;
    t_compute_deriv += elapsed/1000.;
    #endif

    //--------------------------------------------------------------------------
    // STEP 2: Backward pass, compute optimal control law and cost-to-go

    #ifdef TIMESTUFF
    start = std::chrono::system_clock::now();
    #endif

    bool backPassDone = false;
    while (!backPassDone) {
       // update Vx, Vxx, l, L, dV with back_pass
      diverge = 0;
      diverge = backward_pass();

      if (diverge != 0) {
        // cout << "backward pass diverged" << endl;
        dlambda = std::max(dlambda * lambdaFactor, lambdaFactor);
        lambda  = std::max(lambda * dlambda, lambdaMin);
        if (lambda > lambdaMax) break;
        continue;
      }
      backPassDone = true;
    }

    // check for termination due to small gradient
    gnorm = get_gradient_norm(k, us);
    if (gnorm < tolGrad && lambda < 1e-5) {
      #ifdef SHOWPROGRESS
        cout << "\nSUCCESS: gradient norm < tolGrad\n" << endl;
      #endif
      break;
    }

    #ifdef TIMESTUFF
    now = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    // cout << "backward pass took: " << elapsed/1000. << " seconds." << endl;
    t_backward += elapsed/1000.;
    #endif

    //--------------------------------------------------------------------------
    // STEP 3: Forward pass / line-search to find new control sequence, trajectory, cost

    #ifdef TIMESTUFF
    start = std::chrono::system_clock::now();
    #endif

    bool fwdPassDone = false;
    VecOfVecXd xnew(T+1);
    VecOfVecXd unew(T);
    double alpha;

    // TODO parallelize this - split into separate function?
    // Inputs: Alpha, k, x0, dV

    if (backPassDone) { //  serial backtracking line-search
      for (int i=0; i<Alpha.size(); i++) {
        alpha = Alpha(i);
        
        // VecOfVecXd u_plus_feedforward = add_bias_to_u(us, k, alpha);
        VecOfVecXd u_plus_feedforward = us;
        for(unsigned int i=0; i<us.size(); i++) {
          u_plus_feedforward[i] += k[i]*alpha;
        }
        
        new_cost = forward_pass(x0, u_plus_feedforward);
        dcost    = cost_s - new_cost;
        expected = -alpha * (dV(0) + alpha*dV(1));

        if (expected>0) {
          z = dcost/expected;
        }
        else {
          z = sgn(dcost);
          cout << "Warning: non-positive expected reduction. This should not occur" << endl;
        }

        if(z > zMin) {
          fwdPassDone = true;
          break;
        }
      }

      if (!fwdPassDone) {
        // cout << "Forward pass failed" << endl;
        alpha = 0.0; // signals failure of forward pass
      }
    }

    #ifdef TIMESTUFF
    now = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    // cout << "forward pass took: " << elapsed/1000. << " seconds." << endl;
    t_forward += elapsed/1000.;
    #endif

    //--------------------------------------------------------------------------
    // STEP 4: accept step (or not), print status
    #ifdef SHOWPROGRESS
      if (iter==0)
       std::cout << "iteration\tcost\t\treduction\texpect\t\tgrad\t\tlog10(lambda)\n";
    #endif

     if (fwdPassDone) {
      #ifdef SHOWPROGRESS
        printf("%-12d\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                  iter, new_cost, dcost, expected, gnorm, log10(lambda));
      #endif

       // decrease lambda
       dlambda   = min(dlambda / lambdaFactor, 1/lambdaFactor);
       lambda    = lambda * dlambda * (lambda > lambdaMin);

       // accept changes
       cost_s          = new_cost;
       flgChange       = true;

       //terminate?
       if (dcost < tolFun) {
        #ifdef SHOWPROGRESS
           std::cout << "\nSUCCESS: cost change < tolFun\n";
        #endif
         break;
       }
     }
     else { // no cost improvement
      xs = x_old;
      us = u_old;

       // increase lambda
       dlambda  = max(dlambda * lambdaFactor, lambdaFactor);
       lambda   = max(lambda * dlambda, lambdaMin);

      #ifdef SHOWPROGRESS
        printf("%-12d\t%-12s\t%-12.3g\t%-12.3g\t%-12.3g\t%-12.1f\n",
                  iter,"NO STEP", dcost, expected, gnorm, log10(lambda));
      #endif

       // terminate?
       if (lambda > lambdaMax) {
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

  #ifdef TIMESTUFF
    auto all_end = std::chrono::system_clock::now();
    long int t_total = std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_start).count();
    cout << "Time breakdown: " << t_total/1000. << endl;
    cout << "compute_derivatives: " << t_compute_deriv << endl;
    cout << "backward pass: " << t_backward << endl;
    cout << "forward pass: " << t_forward << endl;
    cout << "other stuff: " << t_total/1000. - (t_compute_deriv+t_backward+t_forward) << endl;
  #endif

  output_to_csv("ilqr_result.csv");

} //generate_trajectory


// Saves new state and control sequence in xs, us
//TODO make inputs/outputs explicit here?
double iLQR::forward_pass(const VectorXd &x0, const VecOfVecXd &u) {
  double total_cost = 0;

  VectorXd x_curr = x0;
  VectorXd x1;
  VectorXd u_curr;

  VecOfVecXd x_new(T+1);
  x_new[0] = x0;
  
  for(int t=0; t<T; t++) {
    u_curr = u[t];
    if (K.size()>0) u_curr += K[t]*(x_new[t] - xs[t]); //apply LQR control gains after first iteration

    us[t] = clamp_to_limits(u_curr, model->u_min, model->u_max);
    x1 = model->integrate_dynamics(x_curr, u_curr, dt);
    total_cost += model->cost(x_curr, u_curr);

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
int iLQR::backward_pass() {

  //cost-to-go at end
  Vx[T] = cx[T];
  Vxx[T] = cxx[T];

  dV.setZero();

  for (int i=(T-1); i>=0; i--) { // back up from end of trajectory
    Qx  = cx[i] + (fx[i].transpose() * Vx[i+1]);
    Qu  = cu[i] + (fu[i].transpose() * Vx[i+1]);
    Qxx = cxx[i] + (fx[i].transpose() * Vxx[i+1] * fx[i]);
    Qux = cxu[i].transpose() + (fu[i].transpose() * Vxx[i+1] * fx[i]);
    Quu = cuu[i] + (fu[i].transpose() * Vxx[i+1] * fu[i]);

    // Similar to equations 10a and 10b in [Tassa 2012]. Note that regularization is different
    Qux_reg = cxu[i].transpose() + (fu[i].transpose() * Vxx[i+1] * fx[i]);
    QuuF = cuu[i] + (lambda*MatrixXd::Identity(model->u_dims,model->u_dims)) + (fu[i].transpose() * Vxx[i+1] * fu[i]);

    boxQPResult res = boxQP(QuuF, Qu, k[min(i+1,T-1)], model->u_min - us[i], model->u_max - us[i]);

    if(res.result < 1) return i;

    k_i = res.x_opt;
    VectorXi v_free = res.v_free;

    K_i.setZero();
    if (v_free.any()) {
      MatrixXd R = res.H_free;
      MatrixXd Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, v_free));

      int row_i = 0;
      for(int k=0; k<v_free.size(); k++) {
        if(v_free(k)) K_i.row(k) = Lfree.row(row_i++);
      }
    }

    // Update cost-to-go approximation. Equations 11 in [Tassa 2012]
    dV(0) += k_i.transpose()*Qu;
    dV(1) += 0.5*k_i.transpose()*Quu*k_i;

    Vx[i]  = Qx  + K_i.transpose()*Quu*k_i + K_i.transpose()*Qu + Qux.transpose()*k_i;
    Vxx[i] = Qxx + K_i.transpose()*Quu*K_i + K_i.transpose()*Qux + Qux.transpose()*K_i;
    Vxx[i] = 0.5 * (Vxx[i] + Vxx[i].transpose());

    // save controls/gains
    k[i] = k_i;
    K[i] = K_i;
  }

  return 0;
}

// Replaces this line from matlab: g_norm = mean(max(abs(l) ./ (abs(u)+1), [], 1));
// TODO make this cleaner
double iLQR::get_gradient_norm(const VecOfVecXd& l, const VecOfVecXd& u)
{
  std::vector<double> vals(l.size());
  for (unsigned int i=0; i<l.size(); i++) {
    VectorXd v = l[i].cwiseAbs().array() / (u[i].cwiseAbs().array() + 1);
    vals[i] = v.maxCoeff();;
  }
  return std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
}


void iLQR::output_to_csv(const std::string filename) {
  FILE *XU = fopen(filename.c_str(), "w");

  for(unsigned int i=1; i<=xs[0].size(); i++) fprintf(XU, "x%d, ", i);
  for(unsigned int j=0; j<us[0].size(); j++) fprintf(XU, "u%d, ", j);
  fprintf(XU, "u%d\n", int(us[0].size()));

  for(int t=0; t<T; t++) {
    for(unsigned int i=0; i<xs[t].size(); i++) fprintf(XU, "%f, ", xs[t](i));
    for(unsigned int j=0; j<us[t].size()-1; j++) fprintf(XU, "%f, ", us[t](j));
    fprintf(XU, "%f\n", us[t](us[t].size()-1));
  }

  for(unsigned int i=0; i<xs[T].size(); i++) fprintf(XU, "%f, ", xs[T](i));

  fclose(XU);
  cout << "Saved iLQR result to " << filename << endl;
}
