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

#define BOOST_TEST_MODULE test_depth_utils // NOLINT

#include "libam_support/utils/depth_utils.hh"

#include <boost/test/unit_test.hpp>

using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_depth_utils_1)
{
    auto [npr, nnr] = calculate_num_reads(225, 125, 5.0, 5.0, 1);
    BOOST_TEST(9 == npr);
    BOOST_TEST(9 == nnr);
}

BOOST_AUTO_TEST_CASE(test_depth_utils_2)
{
    auto [npr, nnr] = calculate_num_reads(225, 125, 5.0, 5.0, 2);
    BOOST_TEST(10 == npr);
    BOOST_TEST(8 == nnr);
}
