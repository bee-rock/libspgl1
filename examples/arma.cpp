#include <iostream>
#include <memory>
#include <armadillo>
#include <libspgl1/libspgl1.h>

// Client code
namespace libspgl1 {
namespace matrix {

	template<>
    arma::mat multiply(const arma::mat& a, const arma::mat& b)
    {
        return a * b;
    }

	template<>
    arma::mat add(const arma::mat& a, const arma::mat& b)
    {
        return a + b;
    }

	template<>
    size_t n_rows(const arma::mat& a)
    {
        return a.n_rows;
    }

	template<>
    size_t n_cols(const arma::mat& a)
    {
        return a.n_cols;
    }

} // matrix
} // project
// Client code end

int main()
{
    arma::mat a = arma::eye<arma::mat>(3,3);
    arma::mat b = arma::eye<arma::mat>(3,3);
    libspgl1::do_math(a, b);
}
