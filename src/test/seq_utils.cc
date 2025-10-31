/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#define BOOST_TEST_MODULE test_seq_utils // NOLINT

#include "libam_support/Dtypes.h"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>

using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_seq_utils_1)
{
    const std::string aligned_ref("---AAACCCTTTGG-GAA");
    const std::string aligned_query("TAAAAA-CCCA--GTG--");
    const PairwiseAlignment pwa("", "", "", "", "", std::vector<am_qual_t>(), aligned_query, aligned_ref, 0, true);
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
