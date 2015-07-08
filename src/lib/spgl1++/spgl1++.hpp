#include <algorithm>
#include <limits>

namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType& x0){
	libspgl1::Parameters parameters = libspgl1::Parameters(b);
	VectorType sign_x = libspgl1::vector::sign<VectorType>(x0);
	VectorType x = libspgl1::projectI<VectorType>(x0, parameters.tau);
	VectorType r = libspgl1::initialization::compute_r(A, b, x);
	double f = libspgl1::initialization::compute_f(r);
	VectorType g = libspgl1::initialization::compute_g(At, r);
	VectorType dx = libspgl1::projectI<VectorType>(x-g, parameters.tau)-x;
	double dxNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(dx));
	double gStep = libspgl1::initialization::compute_gstep(dxNorm,
														   parameters.stepmin,
														   parameters.stepmax);
	std::cout << "f: " << f << std::endl;
	std::cout << "gstep: " << gStep << std::endl;
	double fBest     = f;
	VectorType xBest = x;
	double fOld      = f;

	double sigma = 0.0;
	double bNorm = libspgl1::math::norm<double>(b, 2.0);
    size_t nNewton  = 0;
    size_t nPrevVals = 10;
    VectorType lastFv(nPrevVals);
    for (size_t i=0; i<nPrevVals; i++){
    	libspgl1::vector::set_element(lastFv, i, std::numeric_limits<int>::max());
    }

	while(true){
		double gNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(-g));
		std::cout << "gNorm: " << gNorm << std::endl;
		double rNorm = libspgl1::math::norm<double>(r, 2);
		std::cout << "rNorm: " << rNorm << std::endl;
		VectorType residual_minus_measurements = r-b;

		double gap = libspgl1::vector::dot<double>(r, residual_minus_measurements) + parameters.tau*gNorm;
		std::cout << "gap: " << gap << std::endl;

		double rGap = abs(gap) / std::max(1.0,f);
		std::cout << "rGap: " << rGap << std::endl;

		double aError1 = rNorm - sigma;
		std::cout << "aError1: " << aError1 << std::endl;

		double aError2 = f - std::pow(sigma, 2) / 2.0;
		std::cout << "aError2: " << aError2 << std::endl;

		double rError1 = std::abs(aError1) / std::max(1.0,rNorm);
		std::cout << "rError1: " << rError1 << std::endl;

		double rError2 = std::abs(aError2) / std::max(1.0,f);
		std::cout << "rError2: " << rError2 << std::endl;

		if (rGap <= std::max(parameters.optTol, rError2) || rError1 <= parameters.optTol ){
			if (rNorm       <=   parameters.bpTol * bNorm){
				break;
			}
			if (rError1     <=  parameters.optTol){
				break;
			}
			if (rNorm       <=  sigma){
				break;
			}
		}
		bool testRelChange1 = (std::abs(f - fOld) <= parameters.decTol * f);
		bool testRelChange2 = (std::abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
		bool testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) ||
							   (testRelChange2 && rNorm <= 2 * sigma)) && !testUpdateTau;

	    if(testUpdateTau){
	          double tauOld   = parameters.tau;
	          parameters.tau  = std::max(0.0, parameters.tau + (rNorm * aError1) / gNorm);
	          nNewton  = nNewton + 1;
	          if (parameters.tau < tauOld){
	             x = libspgl1::projectI<VectorType>(x, parameters.tau);

	             r = b - A*x;
	             g =   - At*r;
	             f = libspgl1::vector::dot<double>(r,r) / 2.0;
	             libspgl1::vector::set_element(lastFv, 1, f);
	          }
	    }
		break;
	}

    return dx;
}

} // libspgl1
