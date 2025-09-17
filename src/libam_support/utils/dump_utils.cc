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

#include "art_modern_config.h"

#include "libam_support/utils/dump_utils.hh"

#include <boost/filesystem/operations.hpp>

#ifdef WITH_BOOST_STACKTRACE
#include <boost/stacktrace/safe_dump_to.hpp>
#include <boost/stacktrace/stacktrace.hpp>
#endif

#include <csignal>
#include <fstream>
#include <iostream>

namespace labw::art_modern {
#ifdef WITH_BOOST_STACKTRACE

constexpr char DUMP_FILENAME[] = "./backtrace.dump";
namespace {
    // FIXME: The signal handler is not guaranteed to be thread safe. The return values are also unprocessed.
    void my_signal_handler(const int signum)
    {
        std::signal(signum, SIG_DFL);
        boost::stacktrace::safe_dump_to(DUMP_FILENAME);
        std::raise(SIGABRT);
    }

} // namespace

void handle_dumps()
{
    std::signal(SIGSEGV, &my_signal_handler);
    std::signal(SIGABRT, &my_signal_handler);
    if (boost::filesystem::exists(DUMP_FILENAME)) {
        // there is a backtrace
        std::ifstream ifs(DUMP_FILENAME);

        const boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
        std::cout << "Previous run crashed:\n" << st << std::endl;

        // cleaning up
        ifs.close();
        boost::filesystem::remove(DUMP_FILENAME);
    }
}
#else
void handle_dumps() { }
#endif
} // namespace labw::art_modern