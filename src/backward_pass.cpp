#include "ilqr.h"

/*
 	INPUTS
 	   cx: 2x(T+1)					cu: 2x(T+1)
 		 cuu: nxnx(T+1)				cxx: nxnx(T+1)	cuu: 2x2x(T+1)
 		 fx: nxnx(T+1)				fu: nx2x(T+1)		fxx: none
 		 fxu: None						fuu: none				u: 2xT
    OUTPUTS
 	   Vx: nx(T+1)			Vxx: nxnx(T+1)			k:mxT
 	   K: mxnxT 				dV: 2x1
 		 diverge - returns 0 if it doesn't diverge, timestep where it diverges otherwise
*/
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
			MatrixXd Lfree(m,n);
			Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, v_free));

			// TODO fix this hack
			if(v_free[0]==1 && v_free[1]==1)
				K_i = Lfree;
			else if (v_free[0]==1)
				K_i.row(0) = Lfree;
			else if (v_free[1]==1)
				K_i.row(1) = Lfree;
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
