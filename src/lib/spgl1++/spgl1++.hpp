#pragma once
#include <algorithm>
#include <limits>
#include "spgLineCurvy.hpp"
#include "spgLine.hpp"

namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType& x0){
	std::cout.precision(15);
	libspgl1::Parameters parameters = libspgl1::Parameters(b);
	VectorType sign_x = libspgl1::vector::sign<VectorType>(x0);
	VectorType x = libspgl1::projectI<VectorType>(x0, parameters.tau);
	auto xOld = x;
	VectorType r = libspgl1::initialization::compute_r(A, b, x);
	auto rOld = r;
	double f = libspgl1::initialization::compute_f(r);
	VectorType g = libspgl1::initialization::compute_g(At, r);
	auto gOld = g;
	VectorType dx = libspgl1::projectI<VectorType>(x-g, parameters.tau)-x;
	double dxNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(dx));
	double gStep = libspgl1::initialization::compute_gstep(dxNorm,
														   parameters.stepmin,
														   parameters.stepmax);
	double fBest     = f;
	VectorType xBest = x;
	double fOld      = f;
	size_t maxLineErrors = 10;
	double sigma = 0.0;
	double bNorm = libspgl1::math::norm<double>(b, 2.0);
    size_t nNewton  = 0;
    size_t nPrevVals = 3;
    size_t nLineTot = 0;

    VectorType lastFv(nPrevVals);
    for (size_t i=0; i<nPrevVals; i++){
    	libspgl1::vector::set_element<double>(lastFv, i, -std::numeric_limits<double>::max());
    }
    libspgl1::vector::set_element<double>(lastFv, 0.0, f);
    double tauOld{0};
    bool exit = false;
    int iter = 0;
	//while(true){
    for(size_t i=0;i<parameters.outer_iterations;++i){
		double gNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(-g));
		double rNorm = libspgl1::math::norm<double>(r, 2);
		VectorType residual_minus_measurements = r-b;
		double gap = libspgl1::vector::dot<double>(r, residual_minus_measurements) + parameters.tau*gNorm;
		double rGap = std::abs(gap) / std::max(1.0,f);
		double aError1 = rNorm - sigma;
		double aError2 = f - std::pow(sigma, 2) / 2.0;
		double rError1 = std::abs(aError1) / std::max(1.0, rNorm);
		double rError2 = std::abs(aError2) / std::max(1.0, f);
		if (rGap <= std::max(parameters.optTol, rError2) || rError1 <= parameters.optTol ){
			if (rNorm       <=   parameters.bpTol * bNorm){
				exit = true;
			}
			if (rError1     <=  parameters.optTol){
				exit = true;
			}
			if (rNorm       <=  sigma){
				exit = true;
			}
		}
		bool testRelChange1 = (std::abs(f - fOld) <= parameters.decTol * f);
		bool testRelChange2 = (std::abs(f - fOld) <= 0.1 * f * (std::abs(rNorm - sigma)));
		bool testUpdateTau  = ((testRelChange1 && rNorm >  2.0 * sigma) ||
							   (testRelChange2 && rNorm <= 2.0 * sigma)) && !testUpdateTau;

	    if(testUpdateTau){
	          tauOld   = parameters.tau;
	          parameters.tau  = std::max(0.0, parameters.tau + (rNorm * aError1) / gNorm);
	          nNewton  = nNewton + 1;
	          if (parameters.tau < tauOld){
	             x = libspgl1::projectI<VectorType>(x, parameters.tau);
	         	 r = libspgl1::initialization::compute_r(A, b, x);
	         	 f = libspgl1::initialization::compute_f(r);
	         	 g = libspgl1::initialization::compute_g(At, r);
	             for (size_t i=0; i<nPrevVals; i++){
	             	libspgl1::vector::set_element<double>(lastFv, i, -std::numeric_limits<double>::max());
	             }
	             libspgl1::vector::set_element<double>(lastFv, 0.0, f);
	          }
	    }

	    if(exit){
	    	break;
	    }

	    iter = iter + 1;
	    xOld = x;
	    fOld = f;
	    gOld = g;
	    rOld = r;
	    VectorType g_tmp  = gStep*g;
	    auto v = libspgl1::spgLineCurvy<MatrixType, VectorType>(
	    		A, x, parameters, b, g_tmp, libspgl1::math::max<double>(lastFv));

	    f = v.fNew;
	    x = v.xNew;
	    r = v.rNew;

	    nLineTot = nLineTot + v.iter;
	   if(~v.EXIT_CONVERGED){
		  x    = xOld;
		  f    = fOld;
		  r    = rOld;
		  dx   = libspgl1::projectI<VectorType>(x - gStep*g, parameters.tau) - x;
		  double gtd  = libspgl1::vector::dot<double>(g,dx);
		  auto v1 = libspgl1::spgLine<MatrixType, VectorType>(
				 A, x, b, f, libspgl1::math::max<double>(lastFv), gtd, dx);
		    f = v1.fNew;
		    x = v1.xNew;
		    r = v1.rNew;
		  nLineTot = nLineTot + v1.iter;
	   }

       if (v.EXIT_CONVERGED){
    	   g = libspgl1::initialization::compute_g(At, r);
    	   VectorType s    = x - xOld;
    	   VectorType y    = g - gOld;
    	   double sts  = libspgl1::vector::dot<double>(s,s);
    	   double sty  = libspgl1::vector::dot<double>(s,y);
    	   if(sty <= 0){
    		   gStep = parameters.stepmax;
    	   } else{
    		   gStep = std::min(parameters.stepmax,
    				            std::max(parameters.stepmin, sts/sty));
    	   }
       }else{
          gStep = std::min(parameters.stepmax, gStep);
       }

       if(f > (std::pow(sigma,2) / 2.0)){
          lastFv((iter % nPrevVals)) = f;
          if(fBest > f){
             fBest = f;
             xBest = x;
          }
       }
	}

    return x;
}

} // libspgl1
