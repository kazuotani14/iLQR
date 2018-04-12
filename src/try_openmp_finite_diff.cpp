/*
/usr/local/bin/g++-7 -std=c++11 -I ../include/eigen -fopenmp try_openmp_finite_diff.cpp -o try_openmp_finite_diff
*/

#define EIGEN_DONT_PARALLELIZE

#include "omp.h"
#include <chrono>
#include <iostream>
#include <functional>
#include <vector>

#include "eigen/Eigen/Core"

using std::cout;
using std::endl;

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std::chrono;

inline void finite_diff_vecvec2scalar(std::function<double(VectorXd, VectorXd)> f, const VectorXd& x1, const VectorXd& x2, MatrixXd& out) {
  // assume out is already right size
  VectorXd p1, p2, m1, m2;
  double eps = 1e-4;

  // TODO can we assume symmetric and only calculate upper triangle? 
  for (int i=0; i<x1.size(); i++){
    for (int j=i; j<x2.size(); j++){ 
        p1 = m1 = x1;
        p2 = m2 = x2;
        p1(i) += eps;
        m1(i) -= eps;
        p2(j) += eps;
        m2(j) -= eps;
        out(i,j) = out(j,i) = (f(p1, p2) - f(m1, p2) - f(p1, m2) + f(m1, m2)) / (4*eps*eps);
    }
  }
}

int main(int argc, char *argv[]) {

  // https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
  Eigen::initParallel();

  // make a quadratic function of N_DIMS
  // TODO find typical values for our problem
  int N_DIMS = 10;
  int T = 500;

  MatrixXd Q(N_DIMS, N_DIMS);
  Q.setIdentity();
  Q(0,0) = 2;

  // Set up query points
  VectorXd x0(N_DIMS); 
  x0.setZero();
  x0(0) = 5;

  std::vector<VectorXd> xs(T);
  for(int i=0; i<xs.size(); ++i) xs[i] = x0;

  // calculate finite difference of f(x)=x'Hx at all points in xs
  std::function<double(VectorXd)> f = [&Q](VectorXd x){return x.transpose()*Q*x;};
  std::function<double(VectorXd, VectorXd)> f2 = [&Q](VectorXd x1, VectorXd x2){return x2.transpose()*Q*x1;};

  double eps = 1e-4;

  cout << "Starting finite diff Hessian" << endl;
  auto start = std::chrono::system_clock::now();
  std::vector<MatrixXd> Hs(T);
  for(int i=0; i<Hs.size(); ++i) Hs[i] = MatrixXd::Zero(N_DIMS,N_DIMS);

    
  //#define NUM_THREADS 20 // for some reason, 20 works better than 4 (number of cores)
  // omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for // about 2x speedup with this
  for(int t=0; t<xs.size(); ++t) {

    // for (int i=0; i<N_DIMS; i++){
    //   for (int j=i; j<N_DIMS; j++){
    //     VectorXd pp, pm, mp, mm; //plus-plus, plus-minus, ....
    //     pp = pm = mp = mm = xs[t];
    //     pp(i) += eps;
    //     pp(j) += eps;
    //     pm(i) += eps;
    //     pm(j) -= eps;
    //     mp(i) -= eps;
    //     mp(j) += eps;
    //     mm(i) -= eps;
    //     mm(j) -= eps;

    //     Hs[t](i,j) = Hs[t](j,i) = (f(pp) - f(mp) - f(pm) + f(mm)) / (4*eps*eps);
    //   }
    // }

      finite_diff_vecvec2scalar(f2, xs[t], xs[t], Hs[t]); 
  }

  auto now = std::chrono::system_clock::now();
  long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  cout << "Finite differences took: " << elapsed/1000. << " seconds." << endl;

  // for (int i=0; i<Hs.size(); ++i) cout << Hs[i] << "\n---" << endl;
}