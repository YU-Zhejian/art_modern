/**
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
#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/ds/CoverageInfo.hh"

#include <boost/program_options/variables_map.hpp>

#include <cstdlib>
#include <string>
#include <vector>

namespace labw::art_modern {

class ArtIOParams {
public:
    const std::string input_file_name;
    const INPUT_FILE_TYPE art_input_file_type;
    const INPUT_FILE_PARSER art_input_file_parser;
    const CoverageInfo coverage_info;
    const std::size_t parallel;
    const am_readnum_t batch_size;
    const boost::program_options::variables_map vm;
    const std::vector<std::string> args;
};
} // namespace labw::art_modern
