#include "iLQR.h"

MatXd iLQR::rows_w_ind(MatXd &mat, VecXd &indices)
{
	MatXd submat;
  if(mat.rows() != indices.size()){
    std::cout << "mat.rows != indices.size\n";
    return submat;
  }
	for(int i=0; i<indices.size(); i++){
		if(indices(i)>0){
      submat.conservativeResizeLike(MatXd(submat.rows()+1, mat.cols()));
			submat.row(submat.rows()-1) = mat.row(i);
		}
	}
	return submat;
}

int iLQR::backward_pass(const VecOfVecXd &cx, const VecOfVecXd &cu, const VecOfMatXd &cxx, const VecOfMatXd &cxu,
												const VecOfMatXd &cuu, const VecOfMatXd &fx, const VecOfMatXd &fu, const VecOfVecXd &u,

												VecOfVecXd &Vx, VecOfMatXd &Vxx, VecOfVecXd &k, VecOfMatXd &K, Vec2d &dV)
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

	//cost-to-go at end
	Vx[T]     = cx[T];
	Vxx[T]    = cxx[T];

	//initialize Qu, Qx, Qxx, Qxu, Quu, k_i, K_i
	VecXd Qx(n);
	VecXd Qu(m);
	MatXd Qxx(n,n);
	MatXd Qux(m,n);
	MatXd Quu(m,m);
	VecXd k_i(m);
	MatXd K_i(m,n);

	for (int i= T-1; i>=0; i--) // back up from end of trajectory
	{
		Qx  = cx[i]      + fx[i].transpose() * Vx[i+1];
		Qu  = cu[i]      + fu[i].transpose() * Vx[i+1];
		// std::cout << "fu: \n" << fu[i] << '\n';
		// std::cout << "Vx: \n" << Vx[i+1] << '\n';
		// std::cout << "cu: \n" << cu[i] << '\n';

		Qxx = cxx[i] + fx[i].transpose()*Vxx[i+1]*fx[i];

		Qux = cxu[i].transpose() + fu[i].transpose()*Vxx[i+1]*fx[i];

    Quu = cuu[i] + fu[i].transpose()*Vxx[i+1]*fu[i];

		// We are using regularization type 1 (Tassa's code): q_uu+lambda*eye()
		// So Vxx does not have to be regularized
    MatXd Vxx_reg = Vxx[i+1];

		MatXd Qux_reg = cxu[i].transpose() + fu[i].transpose()*Vxx_reg*fx[i];

		// Second term is regularization
    MatXd QuuF = cuu[i] + fu[i].transpose()*Vxx_reg*fu[i] + lambda*Eye2;

		// Impose control limits with boxQP
		VecXd k_i(m);
		MatXd R(m,m);
		VecXd free_v(m);

		// std::cout << "inputs: \n";
		// std::cout << "QuuF: \n" << QuuF << '\n';
		// std::cout << "Qu: \n" << Qu << '\n';
		// std::cout << "k: \n" << k[std::min(i+1,T-1)] << '\n';
		// std::cout << i+1 << ' ' << T-1 << '\n';

		int result = boxQP(QuuF,Qu,k[std::min(i+1,T-1)],  k_i,R,free_v);

		// std::cout << "k_i: \n" << k_i << '\n';

		if (result<1){
			return i;
		}

		if (free_v.any()){
			MatXd Lfree(m,n);
			Lfree = -R.inverse() * (R.transpose().inverse()*rows_w_ind(Qux_reg, free_v));
			// TODO fix this hack
			if(free_v[0]==1 && free_v[1]==1){
				K_i = Lfree;
			}
			else if (free_v[0]==1){
				K_i.row(0) = Lfree;
			}
			else if (free_v[1]==1){
				K_i.row(1) = Lfree;
			}
		}

		// std::cout << "K_i: \n" << K_i << '\n';

		// update cost-to-go approximation
		dV(0) += k_i.transpose()*Qu;
		dV(1) += 0.5*k_i.transpose()*Quu*k_i;

		Vx[i]  = Qx  + K_i.transpose()*Quu*k_i + K_i.transpose()*Qu + Qux.transpose()*k_i;
		Vxx[i] = Qxx + K_i.transpose()*Quu*K_i + K_i.transpose()*Qux + Qux.transpose()*K_i;
		Vxx[i] = 0.5 * (Vxx[i] + Vxx[i].transpose());

		// std::cout << "Vx: \n" << Vx[i] << '\n';
		// std::cout << "Vxx: \n" << Vxx[i] << '\n';
		//
		// std::cout << "\n-------------\n";
		// std::cout << "K_i: \n" << K_i << '\n';
		// getchar();

    // save controls/gains
    k[i]     = k_i;
    K[i]     = K_i;
	}

	return 0;
}
