#pragma once
#include <algorithm>
#include <limits>
#include "spgLineCurvy.hpp"
#include "spgLine.hpp"

namespace libspgl1 {

template<typename MatrixType, typename VectorType>
struct Vectors {
	VectorType x, xBest, xOld, r, rOld, g, gOld, dx, b, lastFv;; //If I declare this first then the code breaks unless I run I compile it in debug mode...what sort of optimization is happening in release mode?

	Vectors(Parameters parameters, const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType &x0) :
		b{b},
		x{libspgl1::projectI<VectorType>(x0, parameters.tau)},
		r{libspgl1::initialization::compute_r(A, b, x)},
		g{libspgl1::initialization::compute_g(At, r)},
		dx{libspgl1::projectI<VectorType>(x-g, parameters.tau)-x},
		xOld{x}, rOld{r}, gOld{g},
		xBest{x},
		lastFv(3)
	{
		for (size_t i=0; i<libspgl1::vector::n_elem(lastFv) ; i++){
			libspgl1::vector::set_element<double>(lastFv, i, -std::numeric_limits<double>::max());
		}
	}
};

template<typename MatrixType, typename VectorType>
struct Quantities {
	double f, dxNorm, gStep, fBest, fOld, sigma, bNorm, tauOld;

	Quantities(const libspgl1::Vectors<MatrixType, VectorType> &vectors, Parameters parameters){
		f = libspgl1::initialization::compute_f(vectors.r);
		dxNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(vectors.dx));
		gStep = libspgl1::initialization::compute_gstep(dxNorm,parameters.stepmin,parameters.stepmax);
		fBest = f;
		fOld  = f;
		sigma = 0.0;
		bNorm = libspgl1::math::norm<double>(vectors.b, 2.0);
		tauOld = 0.0;
	}
};

template<typename MatrixType, typename VectorType>
struct Errors {
	double rNorm, aError1, aError2, rError1, rError2;

	Errors(libspgl1::Vectors<MatrixType, VectorType> vectors, libspgl1::Quantities<MatrixType, VectorType> quantities)
	{
		rNorm = libspgl1::math::norm<double>(vectors.r, 2);
		aError1 = rNorm - quantities.sigma;
		aError2 = quantities.f - std::pow(quantities.sigma, 2) / 2.0;
		rError1 = std::abs(aError1) / std::max(1.0, rNorm);
		rError2 = std::abs(aError2) / std::max(1.0, quantities.f);
	}
};

template<typename MatrixType, typename VectorType>
bool update_tau(libspgl1::Quantities<MatrixType, VectorType> quantities, libspgl1::Parameters parameters, libspgl1::Errors<MatrixType, VectorType> errors){
	bool testRelChange1 = (std::abs(quantities.f - quantities.fOld) <= parameters.decTol * quantities.f);
	bool testRelChange2 = (std::abs(quantities.f - quantities.fOld) <= 0.1 * quantities.f * (std::abs(errors.rNorm - quantities.sigma)));
	return ((testRelChange1 && errors.rNorm >  2.0 * quantities.sigma) || (testRelChange2 && errors.rNorm <= 2.0 * quantities.sigma));
}

template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType &x0)
{
	libspgl1::Parameters parameters(b);
	libspgl1::Vectors <MatrixType, VectorType>vectors(parameters, A, At, b, x0);
	libspgl1::Quantities <MatrixType, VectorType> quantities(vectors, parameters);

	size_t maxLineErrors = 10;
    size_t nNewton  = 0;
    size_t nLineTot = 0;

    libspgl1::vector::set_element<double>(vectors.lastFv, 0.0, quantities.f);
    bool exit = false;
    int iter = 0;
    bool testUpdateTau = true;
	while(true){
		VectorType residual_minus_measurements = vectors.r - b;
		libspgl1::Errors<MatrixType, VectorType> errors(vectors, quantities);
		double gNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(-1.0*vectors.g));
		double gap = libspgl1::vector::dot<double>(vectors.r, residual_minus_measurements) + parameters.tau*gNorm;
		double rGap = std::abs(gap) / std::max(1.0, quantities.f);
		if (rGap <= std::max(parameters.optTol, errors.rError2) || errors.rError1 <= parameters.optTol )
		{
			if (errors.rNorm       <=   parameters.bpTol * quantities.bNorm){ exit = true; }
			if (errors.rError1     <=  parameters.optTol){ exit = true; }
			if (errors.rNorm       <=  quantities.sigma){ exit = true; }
		}

	    if(update_tau(quantities, parameters, errors))
	    {
	    	quantities.tauOld   = parameters.tau;
	    	parameters.tau  = std::max(0.0, parameters.tau + (errors.rNorm * errors.aError1) / gNorm);
	    	nNewton  = nNewton + 1;
	    	if (parameters.tau < quantities.tauOld)
	    	{
	    		vectors.x = libspgl1::projectI<VectorType>(vectors.x, parameters.tau);
	    		vectors.r = libspgl1::initialization::compute_r(A, b, vectors.x);
	    		quantities.f = libspgl1::initialization::compute_f(vectors.r);
	    		vectors.g = libspgl1::initialization::compute_g(At, vectors.r);
	    		for (size_t i=0; i<libspgl1::vector::n_elem(vectors.lastFv); i++)
	    		{
	    			libspgl1::vector::set_element<double>(vectors.lastFv, i, -std::numeric_limits<double>::max());
	    		}
	    		libspgl1::vector::set_element<double>(vectors.lastFv, 0.0, quantities.f);
	    	}
	    }

	    if(exit){
	    	break;
	    }

	    iter = iter + 1;
	    vectors.xOld = vectors.x;
	    quantities.fOld = quantities.f;
	    vectors.gOld = vectors.g;
	    vectors.rOld = vectors.r;
	    VectorType g_tmp  = quantities.gStep*vectors.g;
	    auto projected_line_search = libspgl1::spgLineCurvy<MatrixType, VectorType>(A, vectors.x, parameters, b, g_tmp, libspgl1::math::max<double>(vectors.lastFv));

	    quantities.f = projected_line_search.fNew;
	    vectors.x = projected_line_search.xNew;
	    vectors.r = projected_line_search.rNew;
	    nLineTot = nLineTot + projected_line_search.iter;
	    if(!projected_line_search.EXIT_CONVERGED){
	    	vectors.x    = vectors.xOld;
	    	quantities.f    = quantities.fOld;
	    	vectors.r    = vectors.rOld;
	    	vectors.dx   = libspgl1::projectI<VectorType>(vectors.x - quantities.gStep*vectors.g, parameters.tau) - vectors.x;
	    	double gtd  = libspgl1::vector::dot<double>(vectors.g, vectors.dx);
	    	auto v1 = libspgl1::spgLine<MatrixType, VectorType>(
	    			A, vectors.x, b, quantities.f, libspgl1::math::max<double>(vectors.lastFv), gtd, vectors.dx);
	    	quantities.f = v1.fNew;
	    	vectors.x = v1.xNew;
	    	vectors.r = v1.rNew;
	    	nLineTot = nLineTot + v1.iter;
	    }

	    if (projected_line_search.EXIT_CONVERGED){
	    	vectors.g = libspgl1::initialization::compute_g(At, vectors.r);
	    	VectorType s    = vectors.x - vectors.xOld;
	    	VectorType y    = vectors.g - vectors.gOld;
	    	double sts  = libspgl1::vector::dot<double>(s,s);
	    	double sty  = libspgl1::vector::dot<double>(s,y);
	    	if(sty <= 0){
	    		quantities.gStep = parameters.stepmax;
	    	} else
	    	{
	    		quantities.gStep = std::min(parameters.stepmax, std::max(parameters.stepmin, sts/sty));
	    	}
	    }else
	    {
	    	quantities.gStep = std::min(parameters.stepmax, quantities.gStep);
	    }

	    if(quantities.f > (std::pow(quantities.sigma,2) / 2.0))
	    {
	    	vectors.lastFv((iter % libspgl1::vector::n_elem(vectors.lastFv))) = quantities.f;
	    	if(quantities.fBest > quantities.f)
	    	{
	    		quantities.fBest = quantities.f;
	    		vectors.xBest = vectors.x;
	    	}
	    }
	}

    return vectors.x;
}

} // libspgl1
