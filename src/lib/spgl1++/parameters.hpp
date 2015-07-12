namespace libspgl1 {

struct Parameters {
		size_t max_iterations;
		double bpTol;
		double lsTol;
		double optTol;
		double decTol;
		double stepmin;
		double stepmax;
		bool is_complex;
		double tau;
		int outer_iterations;

		template<typename VectorType>
		Parameters(const VectorType& b) :
			max_iterations(3*libspgl1::vector::n_elem(b)), bpTol(1e-06),
			lsTol(1e-06), optTol(1e-04), decTol(1e-04), stepmin(1e-16),
			stepmax(1e+05), is_complex(false), tau(0.0), outer_iterations(400)
		{}
	};

}
