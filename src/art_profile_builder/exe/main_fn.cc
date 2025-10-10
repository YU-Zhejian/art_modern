#include "art_profile_builder/exe/main_fn.hh"

#include "art_profile_builder/exe/APBConfig.hh"
#include "art_profile_builder/exe/IntermediateEmpDist.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace labw::art_modern {
void view_sam(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
    const std::size_t thread_id, const APBConfig& config)
{
    htsThreadPool tpool = { nullptr, 0 };
    tpool.pool = CExceptionsProxy::assert_not_null(
        hts_tpool_init(config.num_io_threads), "HTSLib", "Failed to init HTS thread pool.");

    auto* in = CExceptionsProxy::assert_not_null(
        hts_open(config.file_path.c_str(), "r"), "HTSLib", "Failed to open HTS file.");
    auto* hdr = CExceptionsProxy::assert_not_null(sam_hdr_read(in), "HTSLib", "Failed to read SAM header.");
    auto* b = CExceptionsProxy::assert_not_null(bam_init1(), "HTSLib", "Failed to init BAM record.");
    am_readnum_t num_valid_reads = 0;
    am_readnum_t num_total_reads = 0;
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &tpool);

    // Skip thread_id reads
    for (std::size_t i = 0; i < thread_id; ++i) {
        if (sam_read1(in, hdr, b) < 0) {
            goto destroy;
        }
    }
    while (true) {
        // read 1 read
        if (sam_read1(in, hdr, b) < 0) {
            goto destroy;
        }
        if ((config.is_pe ? (b->core.flag & BAM_FREAD1 ? ied1 : ied2) : ied1)->parse_read(b)) {
            num_valid_reads++;
        }
        num_total_reads++;
        if (num_total_reads % 500000 == 0) {
            BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                                    << to_si(static_cast<double>(num_total_reads), 2, 1000) << " reads, "
                                    << to_si(static_cast<double>(num_valid_reads), 2, 1000) << " ("
                                    << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads)
                    * 100.0 << ")% valid reads.";
        }
        for (std::size_t i = 1; i < config.num_threads; ++i) {
            // This is really slow...
            if (sam_read1(in, hdr, b) < 0) {
                goto destroy;
            }
        }
    }
destroy:
    BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                            << to_si(static_cast<double>(num_total_reads), 2, 1000) << " reads, "
                            << to_si(static_cast<double>(num_valid_reads), 2, 1000) << " ("
                            << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads) * 100.0
                            << ")% valid reads.";

    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_close(in);
    hts_tpool_destroy(tpool.pool);
}
} // namespace labw::art_modern
