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
#include "libam_support/utils/frontend_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/version_utils.hh"

#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <thread>

namespace po = boost::program_options;

namespace labw::art_modern {
void print_help(const po::options_description& po_desc, const std::string& prefix, const std::string& suffix)
{
    if (!prefix.empty()) {
        std::cout << prefix << std::endl;
    }
    std::cout << po_desc << std::endl;
    if (!suffix.empty()) {
        std::cout << suffix << std::endl;
    }
}

po::variables_map generate_vm_while_handling_help_version(const po::options_description& po_desc, const int argc,
    char** argv, const std::string& prefix, const std::string& suffix)
{
    if (argc == 1) {
        // No command line arguments.
        if (is_on_mpi_main_process_or_nompi()) {
            print_help(po_desc, prefix, suffix);
        }
        abort_mpi();
    }
    po::variables_map vm_;

    try {
        store(parse_command_line(argc, argv, po_desc), vm_);
        notify(vm_);
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        if (is_on_mpi_main_process_or_nompi()) {
            print_help(po_desc, prefix, suffix);
        }
        abort_mpi();
    }

    if (vm_.count(ARG_VERSION) != 0U) {
        if (is_on_mpi_main_process_or_nompi()) {
            print_version();
        }
        exit_mpi();
        std::exit(EXIT_SUCCESS);
    }
    if (vm_.count(ARG_HELP) != 0U) {
        if (is_on_mpi_main_process_or_nompi()) {
            print_help(po_desc, prefix, suffix);
        }
        exit_mpi();
        std::exit(EXIT_SUCCESS);
    }
    return vm_;
}

po::options_description general_options()
{
    po::options_description general_opts("General Options");
    general_opts.add_options()(ARG_HELP, "print out usage information");
    general_opts.add_options()(ARG_VERSION, "display version info");
    return general_opts;
}
/**
 * Extract n_threads from parallel argument.
 */
std::size_t n_threads_from_parallel(const int parallel_arg)
{
    std::size_t parallel = parallel_arg;
    const auto max_threads = std::thread::hardware_concurrency();
    if (parallel_arg == PARALLEL_ALL) {
        parallel = max_threads;
    } else if (parallel_arg == PARALLEL_DISABLE) {
        parallel = 1;
    } else if (parallel_arg < -1) {
        BOOST_LOG_TRIVIAL(fatal) << "parallel (" << parallel << ") must be greater than or equal to -1.";
        abort_mpi();
    } else if (static_cast<std::size_t>(parallel_arg) > max_threads) {
        BOOST_LOG_TRIVIAL(warning) << "parallel (" << parallel
                                   << ") is greater than the "
                                      "maximum number of threads available on the system ("
                                   << max_threads << "). Adjusted to " << max_threads << ".";
        parallel = max_threads;
    }
    return parallel;
}
} // namespace labw::art_modern
