#include "iLQR2.h"

// void iLQR::finite_difference(x, u, matrix, function)
// {
// 	// need function pointer
// 	// or just make separate functions for finite differencing different things, like TGlad did
// }


// This has a weird condition based on number of outputs in MATLAB code
// Split into two functions
double iLQR::get_dynamics_and_cost(const Eigen::VectorXd &x, const Eigen::VectorXd u, Eigen::VectorXd &x1)
{
	// returns cost, modifies dx
		x1 = dynamics(x,u);
		double c = cost(x,u);
		return c;
}

// // //TODO
// double iLQR::compute_derivatives(const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::MatrixXD &fx,
// 												 const Eigen::MatrixXD &fu, const Eigen::Matrix2D &cx, constEigen::Matrix2D &cu,
// 												 const Eigen::MatrixXD &cxx, const Eigen::MatrixXD &cxu, const Eigen::MatrixXD &cuu)
// {
// 	// % state and control indices
// 	// ix = 1:10;
// 	// iu = 11:12;
// 	//
// 	// % dynamics first derivatives
// 	// % J - Jacobian, derivative of states wrt states and control inputs
// 	// % n x (n+m) x T , where n=dim(x), m=dim(u), T=horizon
// 	// xu_dyn  = @(xu) car_dynamics(xu(ix,:),xu(iu,:));
// 	// J       = finite_difference(xu_dyn, [x; u]);
// 	// fx      = J(:,ix,:);
// 	// fu      = J(:,iu,:);
// 	//
// 	// % cost first derivatives
// 	// xu_cost = @(xu) car_cost(xu(ix,:),xu(iu,:));
// 	// J       = squeeze(finite_difference(xu_cost, [x; u]));
// 	// cx      = J(ix,:);
// 	// cu      = J(iu,:);
// 	//
// 	// % cost second derivatives
// 	// xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
// 	// JJ      = finite_difference(xu_Jcst, [x; u]);
// 	// JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
// 	// cxx     = JJ(ix,ix,:);
// 	// cxu     = JJ(ix,iu,:);
// 	// cuu     = JJ(iu,iu,:);
//
// }

// void get_dynamics_derivatives
