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

#include "libam_support/writer/WriterInterface.hh"

#include <htslib/bgzf.h>

#include <string>

namespace labw::art_modern {
class BGZipWriter final : public WriterInterface {
public:
    explicit BGZipWriter(
        const std::string& filename, int level, std::size_t buffer_size, bool write_gzip, std::size_t num_threads);

    ~BGZipWriter() override;
    DELETE_COPY(BGZipWriter)
    DELETE_MOVE(BGZipWriter)

    void write(std::string&& data) override;
    void write(const std::string& data) override;
    void flush() override;
    void close() override;
    [[nodiscard]] std::string name() const override;

private:
    BGZF* gz_file_;
    std::size_t buffer_size_;
};
} // namespace labw::art_modern
