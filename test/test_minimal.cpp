#include <libunittest/all.hpp>
#include <libspgl1/libspgl1.h>
#include <armadillo>
using namespace unittest::assertions;



struct test_minimal : unittest::testcase<> {

    static void run()
    {
        UNITTEST_CLASS(test_minimal)
        UNITTEST_RUN(test_NormL1_primal)
    }

    void test_NormL1_primal()
    {
    	arma::vec v = arma::linspace<arma::vec>(-1, 1, 2);
    	arma::vec ones = arma::ones<arma::vec>(2);
    	assert_equal(2.0, libspgl1::NormL1_primal(v, ones));
    }
};
REGISTER(test_minimal);
