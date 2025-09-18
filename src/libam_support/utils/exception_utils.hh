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

#pragma once

#include "art_modern_config.h"

#ifdef WITH_BOOST_STACKTRACE
#include <boost/exception/all.hpp> // NOLINT
#include <boost/stacktrace/stacktrace.hpp>
#endif

namespace labw::art_modern {
template <class E> [[noreturn]] void throw_with_trace(const E& e)
{
#ifdef WITH_BOOST_STACKTRACE
    throw boost::enable_error_info(e) << boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>(
        boost::stacktrace::stacktrace());
#else
    throw e;
#endif
}
} // namespace labw::art_modern
