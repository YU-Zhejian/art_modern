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

//
// Created by yuzj on 2/22/26.
//

#include "libam_support/writer/BGZipWriter.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <boost/log/trivial.hpp>

#include <htslib/bgzf.h>

#include <string>

namespace labw::art_modern {
BGZipWriter::BGZipWriter(const std::string& filename, const int level, const std::size_t buffer_size,
    const bool write_gzip, const std::size_t num_threads)
    : WriterInterface(filename)
    , buffer_size_(buffer_size)
{
    gz_file_ = bgzf_open(filename.c_str(), (write_gzip ? "wg" : "w" + std::to_string(level)).c_str());
    if (gz_file_ == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open file for writing: " << filename;
        abort_mpi();
    }
    bgzf_set_cache_size(gz_file_, static_cast<int>(buffer_size_));
    bgzf_mt(gz_file_, static_cast<int>(num_threads), 0);
}

void BGZipWriter::close()
{
    if (gz_file_ != nullptr) {
        flush();
        const auto errorno = bgzf_close(gz_file_);
        if (errorno != 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to close BGZF file: " << get_filename() << ", error code: " << errorno;
            abort_mpi();
        }
        gz_file_ = nullptr;
    }
}

std::string BGZipWriter::name() const { return "BGZipWriter"; }

void BGZipWriter::flush()
{
    if (gz_file_ != nullptr) {
        const auto errorno = bgzf_flush(gz_file_);
        if (errorno != 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to flush BGZF file: " << get_filename() << ", error code: " << errorno;
            abort_mpi();
        }
    }
}
void BGZipWriter::write(std::string&& data)
{
    if (gz_file_ != nullptr) {
        auto const data_copy = std::move(data);
        const auto errorno = bgzf_write(gz_file_, data_copy.data(), data_copy.size());
        if (errorno < 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to write to BGZF file: " << get_filename()
                                     << ", error code: " << errorno;
            abort_mpi();
        }
    }
}

void BGZipWriter::write(const std::string& data)
{
    if (gz_file_ != nullptr) {
        const auto errorno = bgzf_write(gz_file_, data.data(), data.size());
        if (errorno < 0) {
            BOOST_LOG_TRIVIAL(error) << "Failed to write to BGZF file: " << get_filename()
                                     << ", error code: " << errorno;
            abort_mpi();
        }
    }
}

BGZipWriter::~BGZipWriter() { BGZipWriter::close(); }
} // namespace labw::art_modern
