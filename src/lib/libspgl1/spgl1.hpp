namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b){


	struct parameters {
		size_t max_iterations;
		double bpTol;
		double lsTol;
		double optTol;
		double decTol;
		double stepmin;
		double stepmax;
		bool is_complex;
		double tau;
		parameters(const VectorType& b) :
			max_iterations(3*libspgl1::vector::n_elem(b)), bpTol(1e-06),
			lsTol(1e-06), optTol(1e-04), decTol(1e-04), stepmin(1e-16),
			stepmax(1e+05), is_complex(false), tau(0)
		{}
	};

	double bNorm  = libspgl1::math::norm(b,2);
	std::cout << bNorm << std::endl;
	//bool subspace = false;
	//double stepG  = 1;

	//VectorType x = project(x, parameters.tau);
	//VectorType r = b - A*x;
	//VectorType g = -At*r;
	//double f     = libspgl1::vector::dot(r,r) / 2;


	//double lastFv    = f;
	//double fBest     = f;
	//VectorType xBest = x;
	//double fOld      = f;
	//return r;


}

} // libspgl1
