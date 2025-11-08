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
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <ceu_check/ceu_check_os_macro.h> // NOLINT: For CEU_ON_POSIX

#include <boost/filesystem/operations.hpp>

#ifdef WITH_BOOST_STACKTRACE
#include <boost/stacktrace/safe_dump_to.hpp>
#include <boost/stacktrace/stacktrace.hpp>
#endif

#ifdef CEU_ON_POSIX
#include <pthread.h>
#include <signal.h>
#endif

#include <csignal>
#include <fstream>
#include <iostream>
#include <string>

namespace labw::art_modern {
#ifdef WITH_BOOST_STACKTRACE
namespace {
    constexpr char DUMP_BASE_FILENAME[] = "./backtrace.dump";

    void my_signal_handler(const int signum) noexcept
    {
        const std::string dump_filename = attach_mpi_rank_to_path(DUMP_BASE_FILENAME, mpi_rank_s());

#ifdef CEU_ON_POSIX
        struct sigaction sa { };
        sa.sa_handler = SIG_DFL;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = 0;
        sigaction(signum, &sa, nullptr);
        boost::stacktrace::safe_dump_to(dump_filename.c_str());
        pthread_kill(pthread_self(), SIGABRT);
#else
        std::signal(signum, SIG_DFL);
        std::raise(SIGABRT);
#endif
    }

} // namespace
#endif
void handle_dumps()
{
#ifdef WITH_BOOST_STACKTRACE
#ifdef CEU_ON_POSIX
    {
        struct sigaction sa { };
        sa.sa_handler = my_signal_handler;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = 0;
        sigaction(SIGSEGV, &sa, nullptr);
        sigaction(SIGABRT, &sa, nullptr);
    }
#else
    std::signal(SIGSEGV, &my_signal_handler);
    std::signal(SIGABRT, &my_signal_handler);
#endif

    std::string const possible_dump_filename = attach_mpi_rank_to_path(DUMP_BASE_FILENAME, "0");
    if (boost::filesystem::exists(possible_dump_filename)) {
        // there is a backtrace
        if (is_on_mpi_main_process_or_nompi()) {
            std::ifstream ifs(possible_dump_filename);
            const boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
            std::cout << "Previous run crashed:\n" << st << std::endl;
            // cleaning up
            ifs.close();
        }
        boost::filesystem::remove(possible_dump_filename);
    }
#endif
}
} // namespace labw::art_modern
