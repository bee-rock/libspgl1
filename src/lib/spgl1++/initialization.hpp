namespace libspgl1 {
namespace initialization {

	template<typename MatrixType, typename VectorType>
	VectorType compute_r(const MatrixType& A,
						 const VectorType& b,
						 const VectorType& x){
		return b - A*x;
	}

	template<typename MatrixType, typename VectorType>
	VectorType compute_g(const MatrixType& At,
						 const VectorType& r){
		return -At*r;
	}

	template<typename VectorType>
	double compute_f(const VectorType& r){
		return libspgl1::vector::dot<double>(r,r) / 2.0;
	}

	double compute_gstep(const double& dxNorm,
						 const double& stepmin,
						 const double& stepmax){
		double gStep;
		if (dxNorm < (1.0 / stepmax)){
		    gStep = stepmax;
		}
		else{
			gStep = std::min(stepmin, std::max(stepmin, 1.0/dxNorm));
		}
		return gStep;
	}

} // initialization
} // libspgl1
