#define BOOST_TEST_MODULE test_seq_utils // NOLINT

#include "libam_support/utils/seq_utils.hh"
#include "libam_support/ds/PairwiseAlignment.hh"

#include <boost/test/unit_test.hpp>

#include <string>

using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_seq_utils_1)
{
    const std::string aligned_ref("---AAACCCTTTGG-GAA");
    const std::string aligned_query("TAAAAA-CCCA--GTG--");
    const PairwiseAlignment pwa("", "", "", "", "", aligned_query, aligned_ref, 0, true);
    BOOST_TEST("3I3=1D2=2X2D1=1I1=2D" == cigar_arr_to_str(pwa.generate_cigar_array(false)));
    BOOST_TEST("3I3M1D4M2D1M1I1M2D" == cigar_arr_to_str(pwa.generate_cigar_array(true)));
}
BOOST_AUTO_TEST_CASE(test_seq_utils_2)
{
    std::string ref_seq("ATCGGCTT");
    revcomp_inplace(ref_seq);
    BOOST_TEST("AAGCCGAT" == ref_seq);
}

BOOST_AUTO_TEST_CASE(test_seq_utils_3)
{
    const std::string ref_seq("ATCGGCTT");
    const auto seq2 = revcomp(ref_seq);
    BOOST_TEST("AAGCCGAT" == seq2);
}
BOOST_AUTO_TEST_CASE(test_seq_utils_4)
{
    std::string ref_seq("ATCGGCTTNNZZ\1");
    normalize_inplace(ref_seq);
    BOOST_TEST("ATCGGCTTNNNNN" == ref_seq);
}
