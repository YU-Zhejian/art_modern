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

#include "libam_support/bam/BamTypes.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/sam.h>

#include <atomic>

namespace labw::art_modern {

class BamLFIO final : public LockFreeIO<bam1_t_uptr> {
public:
    DELETE_MOVE(BamLFIO)
    DELETE_COPY(BamLFIO)
    DEFAULT_DESTRUCTOR(BamLFIO)
    void write(bam1_t_uptr ss) override;
    BamLFIO(std::string name, samFile* fp, const sam_hdr_t* h);

    void flush_and_close() override;

private:
    samFile* fp_;
    const sam_hdr_t* h_;
    std::atomic<bool> closed_ { false };
};

} // namespace labw::art_modern
