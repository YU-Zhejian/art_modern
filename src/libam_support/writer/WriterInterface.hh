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

#include <string>
#include <utility>

namespace labw::art_modern {
class WriterInterface {
public:
    explicit WriterInterface(std::string filename)
        : filename_(std::move(filename))
    {
    }
    virtual ~WriterInterface() = default;
    virtual void write(std::string&& data) = 0;
    virtual void write(const std::string& data) = 0;
    DELETE_COPY(WriterInterface)
    DELETE_MOVE(WriterInterface)
    virtual void close() = 0;
    virtual void flush() = 0;
    [[nodiscard]] std::string get_filename() const { return filename_; }
    [[nodiscard]] virtual std::string name() const = 0;

private:
    std::string filename_;
};
} // namespace labw::art_modern
