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

#include "libam_support/lockfree/SimpleLFIO.hh"

#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/fs_utils.hh"

#include <boost/log/trivial.hpp>

#include <atomic>
#include <fstream>
#include <ios>
#include <memory>
#include <ostream>
#include <string>
#include <utility>

namespace labw::art_modern {

void SimpleLFIO::write(const std::unique_ptr<std::string> value)
{
    if (closed_) {
        return;
    }
    out_ << *value;
    num_bytes_out_ += value->size();
}

SimpleLFIO::SimpleLFIO(std::string name, std::string out_path)
    : LockFreeIO(std::move(name))
    , out_path_(std::move(out_path))
{
    prepare_writer(out_path_);
    BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Writer to '" << out_path_ << "' added.";
    out_ = std::ofstream(out_path_, std::ios::out | std::ios::binary);
}

SimpleLFIO::SimpleLFIO(std::string name, std::string out_path, const std::string& preamble)
    : SimpleLFIO(std::move(name), std::move(out_path))
{
    out_ << preamble;
}

void SimpleLFIO::flush_and_close()
{
    if (closed_) {
        return;
    }
    std::flush(out_);
    out_.close();
    BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Writer to '" << out_path_ << "' closed.";
    closed_ = true;
}
} // namespace labw::art_modern
