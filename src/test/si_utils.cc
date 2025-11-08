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

#define BOOST_TEST_MODULE test_si_utils // NOLINT

#include "libam_support/utils/si_utils.hh"

#include <boost/test/unit_test.hpp>

using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_si_utils_1)
{
    BOOST_TEST("0.00" == to_si(0.0));
    BOOST_TEST("-10.00" == to_si(-10.0));
    BOOST_TEST("1.00K" == to_si(1024.0));
    BOOST_TEST("1.00K" == to_si(1024));
}
