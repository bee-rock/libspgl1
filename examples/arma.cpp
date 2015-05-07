#include <iostream>
#include <memory>
#include <armadillo>
#include <src/libspgl1.h>

// Client code
namespace libspgl1 {
namespace matrix {

    arma::mat multiply(const arma::mat& a, const arma::mat& b)
    {
        return a * b;
    }

    arma::mat add(const arma::mat& a, const arma::mat& b)
    {
        return a + b;
    }

    size_t n_rows(const arma::mat& a)
    {
        return a.n_rows;
    }

    size_t n_cols(const arma::mat& a)
    {
        return a.n_cols;
    }

} // matrix
} // project
// Client code end

int main()
{
    arma::mat a(3, 3);
    arma::mat b(3, 3);
    libspgl1::do_math(a, b);
}
