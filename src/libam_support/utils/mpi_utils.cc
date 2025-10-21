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

#include "art_modern_config.h" // NOLINT: for WITH_MPI, WITH_BOOST_STACKTRACE

#include "libam_support/utils/mpi_utils.hh"

#include "libam_support/utils/log_utils.hh"

#include "libam_support/Constants.hh" // NOLINT: Used in MPI functions

#include <boost/log/trivial.hpp>

#ifdef WITH_BOOST_STACKTRACE
#include <boost/stacktrace/stacktrace.hpp>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <string>

namespace labw::art_modern {

namespace {
    [[noreturn]] [[maybe_unused]] void exception_mpi_is_finalized()
    {
        BOOST_LOG_TRIVIAL(fatal) << "MPI is finalized.";
        abort_mpi(EXIT_FAILURE);
    }
    [[noreturn]] [[maybe_unused]] void exception_mpi_not_available()
    {
        BOOST_LOG_TRIVIAL(fatal) << "MPI is not available.";
        abort_mpi(EXIT_FAILURE);
    }
} // namespace

bool have_mpi() noexcept
{
#ifdef WITH_MPI
    return true;
#else
    return false;
#endif
}

bool is_mpi_finalized()
{
#ifdef WITH_MPI
    int mpi_finalized_flag = 0;
    MPI_Finalized(&mpi_finalized_flag);
    return mpi_finalized_flag != 0;
#else
    exception_mpi_not_available();
#endif
}

void exit_mpi() noexcept
{
    BOOST_LOG_TRIVIAL(info) << "EXIT";
#ifdef WITH_MPI
    if (!is_mpi_finalized()) {
        BOOST_LOG_TRIVIAL(debug) << "Finalizing MPI...";
        MPI_Finalize();
        BOOST_LOG_TRIVIAL(debug) << "MPI finalized.";
    } else {
        BOOST_LOG_TRIVIAL(debug) << "MPI already finalized.";
    }
#endif
    flush_all_sinks();
}

[[noreturn]] void abort_mpi([[maybe_unused]] const int status) noexcept
{
    BOOST_LOG_TRIVIAL(info) << "ABORT";
#ifdef WITH_BOOST_STACKTRACE
    BOOST_LOG_TRIVIAL(info) << "Stacktrace:\n" << boost::stacktrace::stacktrace();
#endif
#ifdef WITH_MPI
    BOOST_LOG_TRIVIAL(debug) << "Sending MPI_ABORT...";
    flush_all_sinks();
    MPI_Abort(MPI_COMM_WORLD, status);
    std::abort();
#else
    BOOST_LOG_TRIVIAL(debug) << "Sending std::abort...";
    flush_all_sinks();
    std::abort();
#endif
}

std::size_t mpi_size()
{
    int size = 1;
#ifdef WITH_MPI
    if (is_mpi_finalized()) {
        exception_mpi_is_finalized();
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return static_cast<std::size_t>(size);
#else
    return size;
#endif
}

void init_mpi([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv)
{
#ifdef WITH_MPI
    MPI_Init(argc, argv);
#endif
}
std::size_t mpi_rank()
{
#ifdef WITH_MPI
    if (is_mpi_finalized()) {
        exception_mpi_is_finalized();
    }
    int rank = MPI_UNAVAILABLE_RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
#else
    exception_mpi_not_available();
#endif
}
std::string mpi_rank_s()
{
#ifdef WITH_MPI
    if (is_mpi_finalized()) {
        return "mpi-finalized";
    }
    int rank = MPI_UNAVAILABLE_RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return std::to_string(rank);
#else
    return "nompi";
#endif
}

bool is_on_mpi_main_process_or_nompi()
{
#ifdef WITH_MPI
    return mpi_rank() == 0;
#else
    return true;
#endif
}
std::string mpi_hostname()
{
#ifdef WITH_MPI
    if (is_mpi_finalized()) {
        return "N/A";
    }
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int name_len = 0;
    MPI_Get_processor_name(hostname, &name_len);
    return std::string(hostname, name_len);
#endif
    return "N/A";
}
} // namespace labw::art_modern
