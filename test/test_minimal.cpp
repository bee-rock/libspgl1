#include <libunittest/all.hpp>
#include <libspgl1/all.hpp>
#include <armadillo>

using namespace unittest::assertions;



struct test_minimal : unittest::testcase<> {

    static void run()
    {
        UNITTEST_CLASS(test_minimal)
        UNITTEST_RUN(test_dot_product)
        UNITTEST_RUN(test_norm)
        UNITTEST_RUN(test_NormL1_primal_with_weighting_equal_to_one)
        UNITTEST_RUN(test_NormL1_primal_with_weights_less_than_one)
        //UNITTEST_RUN(test_NormL1_primal_weights_cannot_be_negative)
        //UNITTEST_RUN(test_NormL1_primal_weights_should_be_less_than_or_equal_to_one)
    }

    void test_dot_product()
    {
    	arma::mat row_vector1(2,1);
    	arma::mat row_vector2(2,1);
    	row_vector1.fill(1.0);
    	row_vector2 = row_vector1;
    	assert_equal(2.0, libspgl1::matrix::dot<float>(row_vector1, row_vector2));
    }


    void test_norm()
    {
    	arma::vec v1 = arma::linspace<arma::vec>(-1, 1, 2);
    	assert_equal(2.0, libspgl1::matrix::norm<float>(v1, 1));
    }

    void test_NormL1_primal_with_weighting_equal_to_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-1, 1, 2);
    	arma::vec weights = arma::ones<arma::vec>(2);
    	assert_equal(2.0, libspgl1::NormL1_primal<float>(v, weights));
    }

    void test_NormL1_primal_with_weights_less_than_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(0.5);
    	assert_equal(2.0, libspgl1::NormL1_primal<float>(v, weights));
    }

    void test_NormL1_primal_weights_cannot_be_negative()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(-0.5);
    	assert_throw<std::invalid_argument>([&v, &weights](){ libspgl1::NormL1_primal<float>(v, weights); }, SPOT);
    }

    void test_NormL1_primal_weights_should_be_less_than_or_equal_to_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(1.5);
    	assert_throw<std::invalid_argument>([&v, &weights](){ libspgl1::NormL1_primal<float>(v, weights); }, SPOT);
    }
};
REGISTER(test_minimal);
