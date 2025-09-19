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

#include "libam_support/CExceptionsProxy.hh"

#include <boost/log/trivial.hpp>
#ifdef WITH_BOOST_STACKTRACE
#include <boost/stacktrace/stacktrace.hpp>
#endif

#include <string>
#include <utility>

namespace labw::art_modern {
const char* CExceptionsProxy::what() const noexcept { return "Error occurred in C library"; }

CExceptionsProxy::CExceptionsProxy(std::string c_lib_name, std::string details)
    : c_lib_name_(std::move(c_lib_name))
    , details_(std::move(details))
{
}
void CExceptionsProxy::log() const
{
    BOOST_LOG_TRIVIAL(fatal) << "Error occurred in C library '" << c_lib_name_ << "' due to '" << details_ << "'";
#ifdef WITH_BOOST_STACKTRACE
    BOOST_LOG_TRIVIAL(fatal) << "Stack trace:";
    BOOST_LOG_TRIVIAL(fatal) << boost::stacktrace::stacktrace();
#endif
}

} // namespace labw::art_modern
