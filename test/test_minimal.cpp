#include <libunittest/all.hpp>
#include <libspgl1/all.hpp>
#include <armadillo>

using namespace unittest::assertions;

namespace libspgl1 {
namespace vector {

    arma::vec abs(const arma::mat& a){
    	return arma::abs(a);
    }

} // matrix
} // project

struct test_minimal : unittest::testcase<> {

    static void run()
    {
        UNITTEST_CLASS(test_minimal)
        UNITTEST_RUN(test_dot_product)
        UNITTEST_RUN(test_norm_l1)
        UNITTEST_RUN(test_norm_l2)
        UNITTEST_RUN(test_NormL1_primal_with_weighting_equal_to_one)
        UNITTEST_RUN(test_NormL1_primal_with_weights_less_than_one)
        UNITTEST_RUN(test_parameters)
        UNITTEST_RUN(test_spgl1_basic_example)
        UNITTEST_RUN(test_projectI)
        //UNITTEST_RUN(test_NormL1_primal_weights_cannot_be_negative)
        //UNITTEST_RUN(test_NormL1_primal_weights_should_be_less_than_or_equal_to_one)
    }

    void test_dot_product()
    {
    	arma::mat row_vector1(2,1);
    	arma::mat row_vector2(2,1);
    	row_vector1.fill(1.0);
    	row_vector2 = row_vector1;
    	assert_equal(2.0, libspgl1::vector::dot<double>(row_vector1, row_vector2));
    }


    void test_norm_l1()
    {
    	arma::vec a = {-1,1};
    	assert_equal(2.0, libspgl1::math::norm<double>(a, 1));
    }

    void test_norm_l2()
    {
    	arma::vec a = {3,4};
    	assert_equal(5.0, libspgl1::math::norm<double>(a, 2));
    }

    void test_NormL1_primal_with_weighting_equal_to_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-1, 1, 2);
    	arma::vec weights = arma::ones<arma::vec>(2);
    	assert_equal(2.0, libspgl1::math::NormL1_primal<double>(v, weights));
    }

    void test_NormL1_primal_with_weights_less_than_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(0.5);
    	assert_equal(2.0, libspgl1::math::NormL1_primal<double>(v, weights));
    }

    void test_NormL1_primal_weights_cannot_be_negative()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(-0.5);
    	assert_throw<std::invalid_argument>([&v, &weights](){ libspgl1::math::NormL1_primal<double>(v, weights); }, SPOT);
    }

    void test_NormL1_primal_weights_should_be_less_than_or_equal_to_one()
    {
    	arma::vec v = arma::linspace<arma::vec>(-2, 2, 2);
    	arma::vec weights = arma::colvec(2);
    	weights.fill(1.5);
    	assert_throw<std::invalid_argument>([&v, &weights](){ libspgl1::math::NormL1_primal<double>(v, weights); }, SPOT);
    }

    void test_parameters()
    {
    	arma::vec a = {3,4};
    	//libspgl1::parameters <arma::vec>parameters(a);
    	//assert_equal(6, parameters.max_iterations);
    }

    void test_spgl1_basic_example(){
    	arma::mat A;
    	A.load("/home/brock/workspace/libspgl1-code/test/A_basic.csv");
    	arma::mat At = A.t();
    	arma::vec x;
    	x.load("/home/brock/workspace/libspgl1-code/test/x0_basic.csv");
    	arma::vec b = A*x;
    	arma::vec x0{libspgl1::matrix::n_cols<arma::mat>(A), arma::fill::zeros};
    	//std::cout << "number of columns" << libspgl1::matrix::n_cols<arma::mat>(A) << std::endl;
    	//std::cout << "x0: " << x0 << std::endl;
    	//std::cout << "number of elements" << libspgl1::vector::n_elem(x0) << std::endl;
    	//auto xnew = libspgl1::spgl1(A,At,b,x0);
    	//std::cout << b << std::endl;
    	//std::cout << libspgl1::math::norm<double>(b, 2.0) << std::endl;
    }

    void test_projectI(){
    	arma::vec x_before_project;
    	arma::vec x_after_project_expected;
    	x_before_project.load("/home/brock/workspace/libspgl1-code/test/x_to_project.csv");
    	x_after_project_expected.load("/home/brock/workspace/libspgl1-code/test/x_after_projection.csv");
    	arma::vec x_after_project = projectI(x_before_project, 100.0, 128);
    	std::cout << "norm1 before" << libspgl1::math::norm<double>(x_before_project, 1.0) << std::endl;
    	std::cout << "norm1 after" << libspgl1::math::norm<double>(x_after_project, 1.0) << std::endl;
    	std::cout << "norm1 expected" << libspgl1::math::norm<double>(x_after_project_expected, 1.0) << std::endl;
    }
};
REGISTER(test_minimal);
