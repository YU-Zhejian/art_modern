#include "PairwiseAlignment.hh"
#include "seq_utils.hh"

#define BOOST_TEST_MODULE test_seq_utils
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_seq_utils_1)
{
    string ref = "---AAACCCTTTGG-GAA";
    string query = "TAAAAA-CCCA--GTG--";
    PairwiseAlignment pwa(query, ref);
    BOOST_TEST("3I3=1D2=2X2D1=1I1=2D" == pwa.generate_cigar(false, false));
    BOOST_TEST("3I3M1D4M2D1M1I1M2D" == pwa.generate_cigar(false, true));
    BOOST_TEST("2D1M1I1M2D4M1D3M3I" == pwa.generate_cigar(true, true));
}
