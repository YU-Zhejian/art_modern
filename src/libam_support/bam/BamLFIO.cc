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

#include "art_modern_config.h"

#include "libam_support/bam/BamLFIO.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/bam/BamTypes.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/hts_utils.hh"

#include <htslib/sam.h>

#include <atomic>
#include <utility>

namespace labw::art_modern {

void BamLFIO::write(const bam1_t_uptr ss)
{
    if (closed_) {
        return;
    }
    CExceptionsProxy::assert_numeric(sam_write1(fp_, h_, ss.get()), USED_HTSLIB_NAME, "Failed to write SAM/BAM record",
        false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
}

BamLFIO::BamLFIO(std::string name, samFile* fp, const sam_hdr_t* h)
    : LockFreeIO(std::move(name))
    , fp_(fp)
    , h_(h)
{
}

void BamLFIO::flush_and_close()
{
    if (closed_) {
        return;
    }
    CExceptionsProxy::assert_numeric(sam_flush(fp_), USED_HTSLIB_NAME, "Failed to flush SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    num_bytes_out_ = hts_tell(fp_);
    CExceptionsProxy::assert_numeric(sam_close(fp_), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    closed_ = true;
}
} // namespace labw::art_modern
