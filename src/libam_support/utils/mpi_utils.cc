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

#include "art_modern_config.h" // NOLINT: for WITH_MPI

#include "libam_support/utils/mpi_utils.hh"

#include "libam_support/Constants.hh" // NOLINT: Used in MPI functions

#include <boost/log/trivial.hpp>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring> // NOLINT: Used in MPI functions
#include <string>

namespace labw::art_modern {
void exit_mpi(const int /*status*/)
{
    BOOST_LOG_TRIVIAL(info) << "EXIT";
#ifdef WITH_MPI
    int mpi_finalized_flag;
    MPI_Finalized(&mpi_finalized_flag);
    if (!mpi_finalized_flag) {
        BOOST_LOG_TRIVIAL(debug) << "Finalizing MPI...";
        MPI_Finalize();
        BOOST_LOG_TRIVIAL(debug) << "MPI finalized.";
    } else {
        BOOST_LOG_TRIVIAL(debug) << "MPI already finalized.";
    }
#endif
}

[[noreturn]] void abort_mpi([[maybe_unused]] const int status)
{
    BOOST_LOG_TRIVIAL(info) << "ABORT";
#ifdef WITH_MPI
    BOOST_LOG_TRIVIAL(debug) << "Sending MPI_ABORT...";
    MPI_Abort(MPI_COMM_WORLD, status);
#endif
    BOOST_LOG_TRIVIAL(debug) << "Sending std::abort...";
    std::abort();
}

void bye_mpi()
{
#ifdef WITH_MPI
    BOOST_LOG_TRIVIAL(info) << "Sending BYE signal...";
    auto bufflen = static_cast<int>(std::strlen(MPI_MESSAGE_BYE));
    void* buffer = std::calloc(bufflen, sizeof(char));
    std::memcpy(buffer, MPI_MESSAGE_BYE, bufflen);
    MPI_Request request;

    MPI_Ibcast(buffer, bufflen, MPI_CHAR, MPI_MAIN_RANK, MPI_COMM_WORLD, &request);
    std::free(buffer);
#endif
}
void init_mpi([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv)
{
#ifdef WITH_MPI
    MPI_Init(argc, argv);
#endif
}
std::string mpi_rank()
{
#ifdef WITH_MPI
    int mpi_finalized_flag;
    MPI_Finalized(&mpi_finalized_flag);
    if (mpi_finalized_flag) {
        return std::to_string(MPI_UNAVAILABLE_RANK);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return std::to_string(rank);
#else
    return "nompi";
#endif
}
} // namespace labw::art_modern
