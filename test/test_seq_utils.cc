#define BOOST_TEST_MODULE test_seq_utils
#include <boost/test/unit_test.hpp>
#include <htslib/sam.h>

#include "PairwiseAlignment.hh"
#include "seq_utils.hh"

using namespace std;
using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_seq_utils_1)
{
    string aligned_ref = "---AAACCCTTTGG-GAA";
    string aligned_query = "TAAAAA-CCCA--GTG--";
    PairwiseAlignment pwa("", "", "", "", "", aligned_query, aligned_ref, 0, true);
    BOOST_TEST("3I3=1D2=2X2D1=1I1=2D" == cigar_arr_to_str(pwa.generate_cigar_array(false)));
    BOOST_TEST("3I3M1D4M2D1M1I1M2D" == cigar_arr_to_str(pwa.generate_cigar_array(true)));
}
