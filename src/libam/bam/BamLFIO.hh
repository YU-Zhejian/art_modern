#pragma once

#include "libam/CExceptionsProxy.hh"
#include "libam/bam/BamTypes.hh"
#include "libam/lockfree/LockFreeIO.hh"
#include "libam/utils/class_macros_utils.hh"

#include <htslib/sam.h>

namespace labw::art_modern {
class BamLFIO : public LockFreeIO<bam1_t_uptr> {
public:
    DELETE_MOVE(BamLFIO)
    DELETE_COPY(BamLFIO)

    void write(bam1_t_uptr ss) override
    {
        CExceptionsProxy::assert_numeric(sam_write1(fp_, h_, ss.get()), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    }
    BamLFIO(samFile* fp, const sam_hdr_t* h)
        : fp_(fp)
        , h_(h)
    {
    }
    ~BamLFIO() override { stop(); };

private:
    samFile* fp_;
    const sam_hdr_t* h_;
};

} // namespace labw::art_modern
