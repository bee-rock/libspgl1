

namespace libspgl1 {
namespace matrix {

template<typename MatrixType>
MatrixType multiply(const MatrixType& a, const MatrixType& b){
	return a * b;
}

template<typename MatrixType>
MatrixType add(const MatrixType& a, const MatrixType& b){
	return a + b;
}

template<typename MatrixType>
size_t n_rows(const MatrixType& a){
	return a.n_rows;
}

template<typename MatrixType>
size_t n_cols(const MatrixType& a){
	return a.n_cols;
}

template<typename MatrixType>
MatrixType dot(const MatrixType& a, const MatrixType& b){
	return a % b;
}

template<typename MatrixType>
float norm(const MatrixType& a, float p){
	float intermediate_value = 0;
	for (auto element : a)
		intermediate_value += pow(abs(element), p);
	return pow(intermediate_value,1/p);
}
}
}

namespace libspgl1 {

template<typename MatrixType>
double NormL1_primal(const MatrixType& a, const MatrixType& b)
{
	int p = 1;
	return libspgl1::matrix::norm(libspgl1::matrix::dot(a, b), p);

}

}
