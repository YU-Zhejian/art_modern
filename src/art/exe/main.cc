/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
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

#include "art_modern_config.h"

#include "art/exe/main_fn.hh"
#include "art/exe/parse_args.hh"

#include "libam_support/utils/dump_utils.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

// Boost timer
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <chrono> // NOLINT
#include <cstdlib>
#include <string>
#include <thread> // NOLINT

using namespace labw::art_modern; // NOLINT

namespace {

void handle_mpi_child()
{
    if (have_mpi()) {
        if (is_on_mpi_main_process_or_nompi()) {
            BOOST_LOG_TRIVIAL(info) << "MPI found! Cross-node parallelism enabled.";
            BOOST_LOG_TRIVIAL(info) << "MPI main process started.";
        } else {
            BOOST_LOG_TRIVIAL(info) << "MPI child process with rank " << mpi_rank() << " started.";
        }
    } else {
        BOOST_LOG_TRIVIAL(info) << "MPI not found! Cross-node parallelism disabled.";
    }
}
} // namespace

int main(int argc, char* argv[])
{
    init_mpi(&argc, &argv);
    // 1st round initialization of a working console logger
    init_logger();
    print_banner();
    init_file_logger();
    handle_mpi_child();
    handle_dumps();

    const auto& [art_params, art_io_params] = parse_args(argc, argv);
    BOOST_LOG_TRIVIAL(info) << "Argument parsing finished. Start generating...";

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    generate_all(art_params, art_io_params);
    BOOST_LOG_TRIVIAL(info) << "All reads generated.";
#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%ws wall, %us user + %ss system = %ts CPU (%p%)");
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
