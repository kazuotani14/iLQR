#include "iLQR2.h"

int iLQR::backward_pass(const VecOfVecXd &cx, const VecOfVecXd &cu, const VecOfMatXd &cxx, const VecOfMatXd &cxu,
												const VecOfMatXd &cuu, const VecOfMatXd &fx, const VecOfMatXd &fu, const VecOfVecXd &u,

												VecOfVecXd Vx, VecOfMatXd Vxx, VecOfVecXd k, VecOfMatXd K, Vec2d dV)
{
// /*
// 	INPUTS
// 	   cx: 2x(T+1)					cu: 2x(T+1)
// 		 cuu: nxnx(T+1)				cxx: nxnx(T+1)	cuu: 2x2x(T+1)
// 		 fx: nxnx(T+1)				fu: nx2x(T+1)		fxx: none
// 		 fxu: None						fuu: none				u: 2xT
//   OUTPUTS
// 	   Vx: nx(T+1)			Vxx: nxnx(T+1)			l:2xT
// 	   L: 2xnxT 				dV: 2x1
// 		 diverge - returns 0 if it doesn't diverge, timestep where it diverges otherwise
//  */

	// preprocessing - resize outputs k, K, Vx, Vxx, dV

	//cost-to-go at end
	//Vx(:,T)     = cx(:,T);
	//Vxx(:,:,T)  = cxx(:,:,T);

	//initialize Qu, Qx, Qxx, Qxu, Quu

	int diverge = 0;
	for (int i= T-1; i>0; i--) // back up from end of trajectory TODO check start index
	{
		// Qu  = cu(:,i)      + fu(:,:,i)'*Vx(:,i+1);
		// Qx  = cx(:,i)      + fx(:,:,i)'*Vx(:,i+1);
		// Qux = cxu(:,:,i)'  + fu(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);

		//  //TODO this stuff


	}



	return diverge;
}
