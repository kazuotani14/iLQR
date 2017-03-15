#include "iLQR2.h"

void iLQR::forward_pass(Eigen::VectorXd &x0, Eigen::Matrix2d &u, Eigen::MatrixXd &L,
                  Eigen::Matrix2d &x, Eigen::Matrix2d &du, Eigen::VectorXd &alpha)
{
/*
	INPUTS
 		x0: nx1			u: 2xT			L: 2xnxT			x: nxT
 		du: 2*T 		alpha: 1x1 or 11x1
	OUTPUTS
		fx: nxnx(T+1)			fu: nx2x(T+1)			fxx: none			fxu: none
		fuu: none					cx: nx(T+1)				cu: 2x(T+1)		cxx: nxnx(T+1)
		cxu: nx2n(T+1)		cuu: 2x2x(T+1)
*/



} //forward_pass
