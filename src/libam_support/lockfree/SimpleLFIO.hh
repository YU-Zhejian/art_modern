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
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <atomic>
#include <fstream>
#include <memory>
#include <string>

namespace labw::art_modern {

class SimpleLFIO : public LockFreeIO<std::unique_ptr<std::string>> {
public:
    DELETE_MOVE(SimpleLFIO)
    DELETE_COPY(SimpleLFIO)
    ~SimpleLFIO() override = default;
    void write(std::unique_ptr<std::string> value) override;

    SimpleLFIO(std::string name, std::string out_path);

    SimpleLFIO(std::string name, std::string out_path, const std::string& preamble);

    void flush_and_close() override;

private:
    std::atomic<bool> closed_ { false };
    std::ofstream out_;
    std::string out_path_;
};
} // namespace labw::art_modern
