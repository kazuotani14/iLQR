#include "ilqr.h"

int iLQR::backward_pass()
{
/*
// 	INPUTS
// 	   cx: 2x(T+1)					cu: 2x(T+1)
// 		 cuu: nxnx(T+1)				cxx: nxnx(T+1)	cuu: 2x2x(T+1)
// 		 fx: nxnx(T+1)				fu: nx2x(T+1)		fxx: none
// 		 fxu: None						fuu: none				u: 2xT
//   OUTPUTS
// 	   Vx: nx(T+1)			Vxx: nxnx(T+1)			l:2xT
// 	   L: 2xnxT 				dV: 2x1
// 		 diverge - returns 0 if it doesn't diverge, timestep where it diverges otherwise
*/

	// TODO
	// % tensor multiplication for DDP terms
	// vectens = @(a,b) permute(sum(bsxfun(@times,a,b),1), [3 2 1]);

	int n = model->x_dims;
	int m = model->u_dims;

	//cost-to-go at end
	Vx[T] = cx[T];
	Vxx[T] = cxx[T];

	//initialize Qu, Qx, Qxx, Qxu, Quu, k_i, K_i
	VectorXd Qx(n), Qu(m);
	MatrixXd Qxx(n,n), Qux(m,n), Quu(m,m);
	VectorXd k_i(m);
	MatrixXd K_i(m,n);

	VecOfVecXd k(T+1); //TODO double check this
	VecOfMatXd K(T+1);
	std::fill(k.begin(), k.end(), VectorXd::Zero(m));
	std::fill(K.begin(), K.end(), MatrixXd::Zero(m,n));

	for (int i=(T-1); i>=0; i--) // back up from end of trajectory
	{
		Qx  = cx[i] + (fx[i].transpose() * Vx[i+1]);
		Qu  = cu[i] + (fu[i].transpose() * Vx[i+1]);

		Qxx = cxx[i] + (fx[i].transpose() * Vxx[i+1] * fx[i]);
		Qux = cxu[i].transpose() + (fu[i].transpose() * Vxx[i+1] * fx[i]);
	    Quu = cuu[i] + (fu[i].transpose() * Vxx[i+1] * fu[i]);

		// We are using regularization type 1 (see Tassa's code): q_uu+lambda*eye()
		// So Vxx does not have to be regularized
	    MatrixXd Vxx_reg = Vxx[i+1];

		MatrixXd Qux_reg = cxu[i].transpose() + (fu[i].transpose() * Vxx_reg * fx[i]);
	    MatrixXd QuuF = cuu[i] + (fu[i].transpose() * Vxx_reg * fu[i]) + (lambda*MatrixXd::Identity(2,2));

		// boxQPResult res = boxQP(QuuF, Qu, k[std::min(i+1,T-1)], model->u_min, model->u_max);
		// int result = res.result;
		// MatrixXd R = res.H_free;
		// VectorXd v_free = res.v_free;

		// if (result < 1){
		// 	return i;
		// }

	// 	if (free_v.any()){
	// 		MatXd Lfree(m,n);
	// 		Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, free_v));
	// 		// TODO fix this hack
	// 		if(free_v[0]==1 && free_v[1]==1){
	// 			K_i = Lfree;
	// 		}
	// 		else if (free_v[0]==1){
	// 			K_i.row(0) = Lfree;
	// 		}
	// 		else if (free_v[1]==1){
	// 			K_i.row(1) = Lfree;
	// 		}
	// 	}

		// update cost-to-go approximation
		// TODO remove assumption of dim==2
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
