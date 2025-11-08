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

BOOST_AUTO_TEST_CASE(test_depth_utils_se_1)
{
    auto [npr, nnr] = calculate_num_reads_se(225, 125, 5.0, 5.0);
    BOOST_TEST(9 == npr);
    BOOST_TEST(9 == nnr);
}

BOOST_AUTO_TEST_CASE(test_depth_utils_se_2)
{
    auto [npr, nnr] = calculate_num_reads_se(225, 125, 10.0, 0.0);
    BOOST_TEST(18 == npr);
    BOOST_TEST(0 == nnr);
}

BOOST_AUTO_TEST_CASE(test_depth_utils_pe_1)
{
    auto [npr, nnr] = calculate_num_reads_pe(225, 125, 125, 5.0, 5.0);
    BOOST_TEST(10 == npr);
    BOOST_TEST(8 == nnr);
}

BOOST_AUTO_TEST_CASE(test_depth_utils_pe_0)
{
    BOOST_TEST(10 == num_base_positive_pe(10, 125, 2, 0));
    BOOST_TEST(125 == num_base_positive_pe(10, 125, 0, 2));

    BOOST_TEST(10 == num_base_negative_pe(10, 125, 0, 2));
    BOOST_TEST(125 == num_base_negative_pe(10, 125, 2, 0));

    BOOST_TEST(135 == num_base_positive_pe(10, 125, 2, 2));
    BOOST_TEST(135 == num_base_negative_pe(10, 125, 2, 2));
}

BOOST_AUTO_TEST_CASE(test_depth_utils_pe_2)
{
    auto [npr, nnr] = calculate_num_reads_pe(225, 10, 140, 0.0, 10.0);
    BOOST_TEST(32 == npr);
    BOOST_TEST(0 == nnr);
}

BOOST_AUTO_TEST_CASE(test_depth_utils_pe_3)
{
    auto [npr, nnr] = calculate_num_reads_pe(225, 10, 140, 10.0, 10.0);
    BOOST_TEST(30 == npr);
    BOOST_TEST(30 == nnr);
    std::cout << num_base_positive_pe(10, 140, npr, nnr) << std::endl;
    std::cout << num_base_negative_pe(10, 140, npr, nnr) << std::endl;
    std::cout << (10.0 + 10.0) * (225) << std::endl;
}
