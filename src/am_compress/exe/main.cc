/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "art_modern_config.h" // NOLINT: For WITH_BOOST_TIMER, WITH_MPI

#include "libam_support/utils/dump_utils.hh"
#include "libam_support/utils/frontend_utils.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/si_utils.hh"
#include "libam_support/writer/WriterDispatcher.hh"

#include <boost/log/trivial.hpp>

#include <boost/filesystem.hpp>

#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <cstdlib>
#include <fstream>
#include <ios>
#include <memory>
#include <vector>

using namespace labw::art_modern;

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);

    if (!is_on_mpi_main_process_or_nompi()) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
    init_logger();
    init_file_logger("am_compress");
    handle_dumps();

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif

    boost::program_options::options_description od {};

    WriterDispatcher::patch_options_for_fmt("compressed", od);
    od.add(general_options());
    od.add_options()("i-file", boost::program_options::value<std::string>(), "Path to the input file.");
    boost::program_options::variables_map vm;
    try {
        vm = generate_vm_while_handling_help_version(od, argc, argv);
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        exit_mpi();
        return EXIT_FAILURE;
    }
    if (vm.empty()) {
        // Help or version requested, already handled.
        exit_mpi();
        return EXIT_SUCCESS;
    }

    const std::string o_file = vm["o-compressed"].as<std::string>();
    prepare_writer(o_file);
    auto writer = WriterDispatcher::parse_args_for_fmt("compressed", vm, o_file);
    std::string in_file = vm["i-file"].as<std::string>();
    if (in_file.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "Input file path is empty.";
        exit_mpi();
        return EXIT_FAILURE;
    }
    if (boost::filesystem::exists(in_file)) {
        if (!boost::filesystem::is_regular_file(in_file)) {
            BOOST_LOG_TRIVIAL(warning) << "Input file path '" << in_file << "' exists but is not a regular file.";
        }
    } else {
        BOOST_LOG_TRIVIAL(fatal) << "Input file '" << in_file << "' does not exist.";
        exit_mpi();
        return EXIT_FAILURE;
    }
    std::fstream input_file(vm["i-file"].as<std::string>(), std::ios::in | std::ios::binary);
    if (!input_file) {
        BOOST_LOG_TRIVIAL(fatal) << "Failed to open input file '" << vm["i-file"].as<std::string>() << "'.";
        exit_mpi();
        return EXIT_FAILURE;
    }
    std::vector<char> buffer(4096);
    while (input_file.read(buffer.data(), buffer.size()) || input_file.gcount() > 0) {
        writer->write(std::string(buffer.data(), input_file.gcount()));
    }
    writer->close();

#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time elapsed: " << t.format();
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
