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
#pragma once

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <cstdlib>
#include <string>

namespace labw::art_modern {
constexpr char ARG_VERSION[] = "version";
constexpr char ARG_HELP[] = "help";
void print_help(
    const boost::program_options::options_description& po_desc, const std::string& prefix, const std::string& suffix);

boost::program_options::variables_map generate_vm_while_handling_help_version(
    const boost::program_options::options_description& po_desc, int argc, char** argv, const std::string& prefix = "",
    const std::string& suffix = "");

boost::program_options::options_description general_options();
/**
 * Extract n_threads from parallel argument.
 */
std::size_t n_threads_from_parallel(int parallel_arg);
} // namespace labw::art_modern
