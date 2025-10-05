/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#define BOOST_TEST_MODULE test_fs_utils // NOLINT

#include "libam_support/utils/fs_utils.hh"

#include <boost/test/unit_test.hpp>

#include <string>

using namespace labw::art_modern;

BOOST_AUTO_TEST_CASE(test_fs_utils_1)
{
    BOOST_TEST(attach_mpi_rank_to_path("/dev/null", "1") == "/dev/null.1");
    BOOST_TEST(attach_mpi_rank_to_path("/dev/.null", "1") == "/dev/.null.1");
    BOOST_TEST(attach_mpi_rank_to_path("/dev/a.null", "1") == "/dev/a.1.null");
    BOOST_TEST(attach_mpi_rank_to_path("/dev/.a.null", "1") == "/dev/.a.1.null");

    BOOST_TEST(attach_mpi_rank_to_path("~/null", "1") == "~/null.1");
    BOOST_TEST(attach_mpi_rank_to_path("~/.null", "1") == "~/.null.1");
    BOOST_TEST(attach_mpi_rank_to_path("~/a.null", "1") == "~/a.1.null");
    BOOST_TEST(attach_mpi_rank_to_path("~/.a.null", "1") == "~/.a.1.null");

    BOOST_TEST(attach_mpi_rank_to_path("./null", "1") == "./null.1");
    BOOST_TEST(attach_mpi_rank_to_path("./.null", "1") == "./.null.1");
    BOOST_TEST(attach_mpi_rank_to_path("./a.null", "1") == "./a.1.null");
    BOOST_TEST(attach_mpi_rank_to_path("./.a.null", "1") == "./.a.1.null");

    BOOST_TEST(attach_mpi_rank_to_path("null", "1") == "null.1");
    BOOST_TEST(attach_mpi_rank_to_path(".null", "1") == ".null.1");
    BOOST_TEST(attach_mpi_rank_to_path("a.null", "1") == "a.1.null");
    BOOST_TEST(attach_mpi_rank_to_path(".a.null", "1") == ".a.1.null");
}
