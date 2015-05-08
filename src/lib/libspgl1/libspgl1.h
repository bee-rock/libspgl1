

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

}
}

namespace libspgl1 {

template<typename MatrixType>
void do_math(const MatrixType& a, const MatrixType& b)
{
    auto c = libspgl1::matrix::multiply(a, b);
    std::cout<<libspgl1::matrix::n_rows(c)<<std::endl;
    std::cout<<libspgl1::matrix::n_cols(c)<<std::endl;
    auto d = libspgl1::matrix::add(a, b);
    std::cout<<libspgl1::matrix::n_rows(d)<<std::endl;
    std::cout<<libspgl1::matrix::n_cols(d)<<std::endl;
    std::cout << c << std::endl;
    std::cout << d << std::endl;
}

}
