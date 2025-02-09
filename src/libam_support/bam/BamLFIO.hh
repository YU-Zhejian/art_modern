#pragma once

#include "art_modern_config.h"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/bam/BamTypes.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/hts_utils.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <atomic>
#include <utility>

namespace labw::art_modern {

class BamLFIO : public LockFreeIO<bam1_t_uptr> {
public:
    DELETE_MOVE(BamLFIO)
    DELETE_COPY(BamLFIO)
    ~BamLFIO() override = default;
    void write(bam1_t_uptr ss) override
    {
        if (closed_) {
            return;
        }
        CExceptionsProxy::assert_numeric(sam_write1(fp_, h_, ss.get()), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    }
    BamLFIO(std::string name, samFile* fp, const sam_hdr_t* h)
        : LockFreeIO<bam1_t_uptr>(std::move(name))
        , fp_(fp)
        , h_(h)
    {
    }
    void flush_and_close() override
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

private:
    samFile* fp_;
    const sam_hdr_t* h_;
    std::atomic<bool> closed_ { false };
};

} // namespace labw::art_modern
