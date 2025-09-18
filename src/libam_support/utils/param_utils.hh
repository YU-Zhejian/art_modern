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

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>
#include <boost/program_options/variables_map.hpp>

#include <exception>
#include <string>

namespace labw::art_modern {

template <typename T> T get_param(const boost::program_options::variables_map& vm, const std::string& name);

} // namespace labw::art_modern

template <typename T>
T labw::art_modern::get_param(const boost::program_options::variables_map& vm, const std::string& name)
{
    try {
        return vm[name].as<T>();
    } catch (const std::exception& exp) {
        BOOST_LOG_TRIVIAL(fatal) << exp.what();
        abort_mpi();
    }
}
