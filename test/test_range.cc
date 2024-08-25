#include "seq_utils.hh"

#define BOOST_TEST_MODULE test_range
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test_range_1)
{
    std::vector<int> vec1;
    std::vector<int> vec2;
    vec1 = labw::art_modern::range(0, 5, 1);
    vec2 = std::vector<int> { 0, 1, 2, 3, 4 };
    BOOST_CHECK_EQUAL_COLLECTIONS(vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
    vec1 = labw::art_modern::range(5, 0, -1);
    vec2 = std::vector<int> { 5, 4, 3, 2, 1 };
    BOOST_CHECK_EQUAL_COLLECTIONS(vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
}
