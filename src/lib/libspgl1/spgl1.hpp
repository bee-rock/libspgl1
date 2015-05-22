namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, VectorType& x0){

	libspgl1::Parameters parameters = libspgl1::Parameters(b);
	VectorType x = libspgl1::projectI<VectorType>(x0, parameters.tau, libspgl1::vector::n_elem(b));
	//VectorType r = libspgl1::initialization::compute_r<VectorType>(A, b, x);
	//double f     = libspgl1::initialization::compute_f(At, r);
	//VectorType g = libspgl1::initialization::compute_g<VectorType>(At, r);
	return x;
}

} // libspgl1
