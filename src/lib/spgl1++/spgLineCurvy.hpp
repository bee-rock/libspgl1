#pragma once

namespace libspgl1{

template<typename VectorType>
struct spgLineCurvyVars {
	double fNew;
	VectorType xNew;
	VectorType rNew;
	double iter;
	double step;
	double err;
	bool EXIT_CONVERGED;

	spgLineCurvyVars(const VectorType &x) :
		fNew(0.0),
		xNew(x),
		rNew(x),
		iter(0),
		step(1.0),
		err(0),
		EXIT_CONVERGED(false)
		{}
};

template <typename MatrixType, typename VectorType>
spgLineCurvyVars<VectorType> spgLineCurvy(const MatrixType &A, VectorType &x, libspgl1::Parameters parameters, const VectorType &b, const VectorType &g, double fMax){

	spgLineCurvyVars <VectorType>v(x);
	double gamma  = 1e-4;
	size_t maxIts = 10;
	double sNorm  =  0.0;
	double gNorm = 0.0;
	double sNormOld = 0.0;
	double scale  =  1.0;      // Safeguard scaling.  (See below.)
	int nSafe  =  0;      // No. of safeguarding steps.
	size_t n = libspgl1::vector::n_elem(x);
	auto s = x;

	while(true){
		v.xNew = libspgl1::projectI<VectorType>(x - v.step*scale*g, parameters.tau);
		v.rNew = libspgl1::initialization::compute_r(A, b, v.xNew);
	    v.fNew = libspgl1::initialization::compute_f(v.rNew);
	    s = v.xNew - x;
	    double gts = scale * libspgl1::vector::dot<double>(g, s);
	    if(gts >= 0){
			std::cout << "no descent" << std::endl;
	    	v.EXIT_CONVERGED = false;
	    	return v;
	    }

	    if(v.fNew < fMax + gamma*v.step*gts){
	    	v.EXIT_CONVERGED = true;
	    	return v;
	    }else if(v.iter >= maxIts){
	    	v.EXIT_CONVERGED = false;
	    	return v;
	    }
	    v.iter = v.iter + 1;
	    v.step = v.step / 2.0;

	    sNormOld  = sNorm;
	    sNorm     = libspgl1::math::norm<double>(s, 2.0) / std::sqrt(static_cast<double>(n));
	    if(std::abs(sNorm - sNormOld) <= 1e-6 * sNorm){
	    	gNorm = libspgl1::math::norm<double>(g, 2.0) / std::sqrt(static_cast<double>(n));
	    	scale = sNorm / gNorm / std::pow(2.0, nSafe);
	        nSafe = nSafe + 1;
	    }
	}
}
}
