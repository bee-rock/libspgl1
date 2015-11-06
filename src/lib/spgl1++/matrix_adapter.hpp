#pragma once
namespace libspgl1 {
namespace matrix {

template<typename MatrixType, typename VectorType>
VectorType matvec(const MatrixType& A, const VectorType& b){
	return A * b;
}

template<typename MatrixType>
size_t n_cols(const MatrixType& a){
	return a.n_cols;
}

template<typename MatrixType>
size_t n_rows(const MatrixType& a){
	return a.n_rows;
}
} // matrix
} // libspgl1
