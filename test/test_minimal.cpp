#include <libunittest/all.hpp>
using namespace unittest::assertions;

struct test_minimal : unittest::testcase<> {

    static void run()
    {
        UNITTEST_CLASS(test_minimal)
        UNITTEST_RUN(test_value_is_true)
    }

    void test_value_is_true()
    {
    	assert_true(true);
    }
};
REGISTER(test_minimal);
