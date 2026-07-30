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

#include <fstream>
#include <string>
#include <vector>

namespace labw::art_modern {
class SimpleWriter : public WriterInterface {
public:
    DELETE_COPY(SimpleWriter)
    DELETE_MOVE(SimpleWriter)
    explicit SimpleWriter(const std::string& filename, std::size_t buffer_size);
    ~SimpleWriter() override;

    void write(std::string&& data) override;
    void close() override;
    void flush() override;
    void write(const std::string& data) override;
    [[nodiscard]] std::string name() const override;

private:
    std::ofstream ofs_;
    std::vector<char> buffer_;
};
} // namespace labw::art_modern
