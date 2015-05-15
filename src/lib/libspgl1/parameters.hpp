template<typename VectorType>
struct parameters {

	size_t iterations = 5000000; //todo
	size_t nPrevVals = 3;
	double bpTol = 1e-06;
	double lsTol = 1e-06;
	double optTol = 1e-04;
	double decTol = 1e-04;
	double stepMin = 1e-16;
	double stepmax = 1e+05;
};

