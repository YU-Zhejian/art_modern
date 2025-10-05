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

#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "libam_support/utils/fs_utils.hh"

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/log/trivial.hpp>

#include <string>

namespace labw::art_modern {
void validate_input_filename(const std::string& input_file_path, const std::string& arg_name)
{
    if (input_file_path.empty()) {
        BOOST_LOG_TRIVIAL(fatal) << "An input file path for --" << arg_name << " must be specified.";
        abort_mpi();
    }
    if (!boost::filesystem::exists(input_file_path)) {
        BOOST_LOG_TRIVIAL(fatal) << "Input file for --" << arg_name << " at '" << input_file_path
                                 << "' does not exist.";
        abort_mpi();
    }
    if (!boost::filesystem::is_regular_file(input_file_path)) {
        BOOST_LOG_TRIVIAL(warning) << "Input file for --" << arg_name << " at '" << input_file_path
                                   << "' is not a regular file.";
#ifdef WITH_MPI
        BOOST_LOG_TRIVIAL(fatal) << "Irregular file is NOT allowed under MPI.";
        abort_mpi();
#endif
    }
}

am_filelen_t get_file_size(const std::string& file_path) noexcept
{
    if (!boost::filesystem::is_regular_file(file_path)) {
        return -1;
    }

    try {
        return static_cast<am_filelen_t>(boost::filesystem::file_size(file_path));
    } catch (const boost::filesystem::filesystem_error&) {
        return -1;
    }
}
void ensure_directory_exists(const std::string& dir_path)
{
    if (boost::filesystem::exists(dir_path)) {
        if (boost::filesystem::is_directory(dir_path)) {
            return;
        }
        BOOST_LOG_TRIVIAL(fatal) << "Path '" << dir_path << "' exists but is not a directory.";
        abort_mpi();
    }
    BOOST_LOG_TRIVIAL(info) << "Boost::filesystem: mkdir -p '" << dir_path << "'";
    boost::filesystem::create_directories(dir_path);
}
void prepare_writer(const std::string& output_file_path)
{
    const auto parent_path = boost::filesystem::path(output_file_path).parent_path().string();
    if (!parent_path.empty()) {
        ensure_directory_exists(parent_path);
    }
    if (boost::filesystem::exists(output_file_path)) {
        if (!boost::filesystem::is_regular_file(output_file_path)) {
#ifdef WITH_MPI
            BOOST_LOG_TRIVIAL(fatal) << "Irregular file is NOT allowed under MPI.";
            abort_mpi();
#else
            BOOST_LOG_TRIVIAL(warning) << "Output file at '" << output_file_path
                                       << "' exists but is not a regular file. Will be overwritten.";
#endif
        } else {
            BOOST_LOG_TRIVIAL(warning) << "Output file at '" << output_file_path << "' exists and will be overwritten.";
        }
    }
}

std::string attach_mpi_rank_to_path(const std::string& file_path, [[maybe_unused]] const std::string& rank_str)
{
#ifdef WITH_MPI
    // Firstly, extract the file name and path
    const auto file_path_path = boost::filesystem::path(file_path);
    const auto filename = file_path_path.filename().string();
    const auto parent_path = file_path_path.parent_path();

    // Find the last dot in the filename
    const auto last_dot_pos = filename.find_last_of('.');

    std::string new_filename;
    if (last_dot_pos == std::string::npos || last_dot_pos == 0) {
        // No dot found or dot is the first character, append rank at the end
        new_filename = filename + "." + rank_str;
    } else {
        // Insert rank before the last dot
        new_filename = filename.substr(0, last_dot_pos) + "." + rank_str + filename.substr(last_dot_pos);
    }

    // Reconstruct the full path
    if (parent_path.empty()) {
        return new_filename; // No parent path, return just the new filename
    } else {
        return (parent_path / new_filename).string(); // Combine parent path and new filename
    }
#else
    return file_path;
#endif
}

} // namespace labw::art_modern
// labw
