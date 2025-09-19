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

#include "profile_utils.hh"

#include <ceu_check/ceu_check_os_macro.h>

#ifdef CEU_ON_POSIX
#include <sys/resource.h>
#endif

#include <boost/log/trivial.hpp>

namespace labw::art_modern::details {
#ifdef CEU_ON_POSIX
void print_memory_usage(const char* file, const int line)
{
    // NOLINTBEGIN
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        BOOST_LOG_TRIVIAL(info) << file << ":" << line << ": " << "Memory usage: " << usage.ru_maxrss * 1024;
    }
    // NOLINTEND
}
#else
void print_memory_usage(const char* /*file*/, const int /*line*/)
{
    // Do nothing on non-POSIX systems
}
#endif

} // namespace lab::art_modern::details
