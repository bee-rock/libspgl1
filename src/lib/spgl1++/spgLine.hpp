#pragma once

namespace libspgl1{

template <typename MatrixType, typename VectorType>
spgLineCurvyVars<VectorType> spgLine(const MatrixType &A, VectorType &x, const VectorType &b, double f, double fMax, double gtd, VectorType &dx){
	spgLineCurvyVars <VectorType>v(x);
	double gamma  = 1e-4;
	size_t maxIts = 10;
	gtd    = -std::abs(gtd);

	while(true){
	    v.xNew = x + v.step*dx;
		v.rNew = libspgl1::initialization::compute_r(A, b, v.xNew);
	    v.fNew = libspgl1::initialization::compute_f(v.rNew);

	    if(v.fNew < fMax + gamma*v.step*gtd){
	    	v.EXIT_CONVERGED = true;
	    	return v;
	    }else if(v.iter >= maxIts){
	    	v.EXIT_CONVERGED = false;
	    	return v;
	    }

	    v.iter = v.iter + 1;

	    if(v.step <= 0.1){
	       v.step  = v.step / 2.0;
	    } else {
	       auto tmp = (-gtd*std::pow(v.step, 2.0)) / (2.0*(v.fNew-f-v.step*gtd));
	       if(tmp < 0.1 || tmp > 0.9*v.step){
	          tmp = v.step / 2.0;
	       }
	       v.step = tmp;
	    }
	}
}
}
