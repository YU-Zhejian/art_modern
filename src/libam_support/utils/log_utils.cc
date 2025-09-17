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

#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "libam_support/utils/log_utils.hh"

#ifdef WITH_MPI
#include "libam_support/utils/mpi_log_attributes.hh"
#endif

#include <boost/filesystem/operations.hpp>
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/current_thread_id.hpp>
#include <boost/log/core/core.hpp>
#include <boost/log/keywords/file_name.hpp>
#include <boost/log/keywords/filter.hpp>
#include <boost/log/keywords/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>
#include <string>

namespace logging = boost::log;
namespace expr = logging::expressions; // NOLINT

namespace labw::art_modern {

void init_logger()
{
    auto core = logging::core::get();
    core->remove_all_sinks();
    core->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    core->add_global_attribute("ThreadID", boost::log::attributes::current_thread_id());
#ifdef WITH_MPI
    core->add_global_attribute("MPIRank", MPIRankLoggerAttribute());
    core->add_global_attribute("MPIHostName", MPIHostNameLoggerAttribute());
#endif
    const auto sink = add_console_log(std::cerr, boost::log::keywords::format = "[%TimeStamp%] %Severity%: %Message%",
        logging::keywords::filter = logging::trivial::severity >= logging::trivial::info);

    core->add_sink(sink);
}
void init_file_logger(const std::string& log_dir)
{
    if (!boost::filesystem::exists(log_dir)) {
        boost::filesystem::create_directories(log_dir);
    }
#ifdef WITH_MPI
    int mpi_num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
    for (int i = 0; i < mpi_num_procs; i++) {
        std::stringstream file_name_ss;
        file_name_ss << log_dir << "/" << i << ".log";
        logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
            logging::keywords::format
            = "[%TimeStamp%] [T=%ThreadID%@MPI=%MPIRank%:%MPIHostName%] %Severity%: %Message%",
            logging::keywords::filter = expr::attr<int>("MPIRank") == i);
    }
    std::stringstream file_name_ss;
    file_name_ss << log_dir << "/" << "nompi" << ".log";
    logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
        logging::keywords::format = "[%TimeStamp%] [T=%ThreadID%@MPI=%MPIRank%:%MPIHostName%] %Severity%: %Message%",
        logging::keywords::filter = expr::attr<int>("MPIRank") == MPI_UNAVAILABLE_RANK);
#else
    std::stringstream file_name_ss;
    file_name_ss << log_dir << "/" << "nompi" << ".log";
    logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
        logging::keywords::format = "[%TimeStamp%] [T=%ThreadID%] %Severity%: %Message%");
#endif
}
} // namespace labw::art_modern