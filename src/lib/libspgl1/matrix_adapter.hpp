namespace libspgl1 {
namespace matrix {

template<typename MatrixType, typename VectorType>
VectorType matvec(const MatrixType& a, const VectorType& b){
	return a * b;
}

} // matrix
} // libspgl1
