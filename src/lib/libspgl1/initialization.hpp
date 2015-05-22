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

	template<typename MatrixType, typename VectorType>
	double compute_f(const MatrixType& At,
					 const VectorType& r){
		return libspgl1::vector::dot<double>(r,r) / 2.0;
	}

} // initialization
} // libspgl1
