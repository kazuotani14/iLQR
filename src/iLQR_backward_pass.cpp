#include "iLQR2.h"

int iLQR::back_pass(Eigen::Matrix2d &cx, Eigen::Matrix2d &cu, Eigen::MatrixXd &cxx,
										Eigen::MatrixXd &cxu, Eigen::MatrixXd &cuu, Eigen::MatrixXd &fx,
										Eigen::MatrixXd &fu, Eigen::Matrix2d &u,

										Eigen::Matrix2d &Vx, Eigen::MatrixXd &Vxx, Eigen::Matrix2d &l,
										Eigen::MatrixXd &L, Eigen::Vector2d &dV)
{
/*
	INPUTS
	   cx: 2x(T+1)					cu: 2x(T+1)
		 cuu: nxnx(T+1)				cxx: nxnx(T+1)	cuu: 2x2x(T+1)
		 fx: nxnx(T+1)				fu: nx2x(T+1)		fxx: none
		 fxu: None						fuu: none				u: 2xT
  OUTPUTS
	   Vx: nx(T+1)			Vxx: nxnx(T+1)			l:2xT
	   L: 2xnxT 				dV: 2x1
		 diverge - returns 0 if it doesn't diverge, timestep where it diverges otherwise
 */

 int diverge = 0;

 //TODO this stuff

 return diverge;
} //backward_pass
