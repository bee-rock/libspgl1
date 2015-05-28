namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType& x0){

	libspgl1::Parameters parameters = libspgl1::Parameters(b);
	VectorType x = libspgl1::projectI<VectorType>(x0, parameters.tau);
	VectorType r = libspgl1::initialization::compute_r(A, b, x);
	double f = libspgl1::initialization::compute_f(r);
	VectorType g = libspgl1::initialization::compute_g(At, r);
	VectorType dx = libspgl1::projectI<VectorType>(x-g, parameters.tau)-x;
	return x0;
}

} // libspgl1
