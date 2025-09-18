/**
 * @brief Argument parser
 *
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

#pragma once
#include "art/lib/ArtIOParams.hh"
#include "art/lib/ArtParams.hh"

#include <tuple>

namespace labw::art_modern {

/**
 * @brief Parse command line arguments and return the parsed parameters.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return Parsed parameters.
 */
std::tuple<ArtParams, ArtIOParams> parse_args(int argc, char** argv);

} // namespace labw::art_modern
