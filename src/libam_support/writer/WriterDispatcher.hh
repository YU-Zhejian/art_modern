/**
 * Copyright 2024-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/compression_utils.hh"

#include "libam_support/writer/WriterInterface.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>

namespace labw::art_modern {

class WriterDispatcher {
public:
    DELETE_COPY(WriterDispatcher)
    DELETE_MOVE(WriterDispatcher)
    DEFAULT_DESTRUCTOR(WriterDispatcher)

    static std::unique_ptr<WriterInterface> writer_dispatch(const std::string& filename,
        CompressionType compression_type = CompressionType::GZIP, int compression_level = 6,
        std::size_t buffer_size = 1 << 20, std::size_t num_threads = 1);

    static std::unique_ptr<WriterInterface> parse_args_for_fmt(
        const std::string& fmt_name, const boost::program_options::variables_map& vm, const std::string& resolved_path);
    static void patch_options_for_fmt(const std::string& fmt_name, boost::program_options::options_description& desc);
};
} // namespace labw::art_modern
