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

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <string>

#include <htslib/faidx.h>
#include <htslib/hts.h>

namespace labw::art_modern {

/**
 * Please note that this method have no lock.
 *
 * That is, all threads MUST own their own FaidxFetch.
 */
class FaidxFetch final : public BaseFastaFetch {
public:
    DELETE_COPY(FaidxFetch)
    DELETE_MOVE(FaidxFetch)

    explicit FaidxFetch(faidx_t* faidx);
    explicit FaidxFetch(const std::string& file_name);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    ~FaidxFetch() override;

private:
    faidx_t* faidx_;
};
} // namespace labw::art_modern
