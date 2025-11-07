#include "art_profile_builder/exe/main_fn.hh"

#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <art_modern_config.h>
#include <cstdlib>
#include <memory>
#include <string>

#define READ_1                                                                                                         \
    retv = sam_read1(in, hdr, b);                                                                                      \
    num_parsed_reads++;                                                                                                \
    if (retv == -1 /** EOF **/) {                                                                                      \
        goto destroy;                                                                                                  \
    } else if (retv < -1) {                                                                                            \
        goto except;                                                                                                   \
    }
constexpr std::size_t REPORT_SIZE = 10000000; // 10 Million reads

namespace labw::art_modern {

void view_sam(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
    const std::size_t thread_id, const APBConfig& config)
{
    int retv = 0;

    htsThreadPool tpool = { nullptr, 0 };
    tpool.pool = CExceptionsProxy::assert_not_null(
        hts_tpool_init(static_cast<int>(config.num_io_threads)), USED_HTSLIB_NAME, "Failed to init HTS thread pool.");

    auto* in = CExceptionsProxy::assert_not_null(
        hts_open(config.input_file_path.c_str(), "r"), USED_HTSLIB_NAME, "Failed to open HTS file.");
    auto* hdr = CExceptionsProxy::assert_not_null(sam_hdr_read(in), USED_HTSLIB_NAME, "Failed to read SAM header.");
    auto* b = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    am_readnum_t num_valid_reads = 0;
    am_readnum_t num_total_reads = 0;
    am_readnum_t num_parsed_reads = 0;
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &tpool);

    // Skip thread_id reads
    for (std::size_t i = 0; i < thread_id; ++i) {
        READ_1
    }
    while (true) {
        // read 1 read
        READ_1
        if (((config.is_pe && ((b->core.flag & BAM_FREAD2) != 0)) ? ied2 : ied1)->parse_read(b)) {
            num_valid_reads++;
        }
        num_total_reads++;
        if (num_total_reads % REPORT_SIZE == 0) {
            BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                                    << to_si(num_total_reads, 2, static_cast<decltype(num_total_reads)>(1000))
                                    << " reads, "
                                    << to_si(num_valid_reads, 2, static_cast<decltype(num_total_reads)>(1000)) << " ("
                                    << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads)
                    * 100.0 << ")% valid reads.";
        }
        // Skip num_threads - 1 reads
        for (std::size_t i = 1; i < config.num_threads; ++i) {
            READ_1
        }
    }
except:
    BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": EXCEPT when parsing read " << num_parsed_reads << ".";
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_close(in);
    hts_tpool_destroy(tpool.pool);
    abort_mpi();
destroy:
    BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                            << to_si(num_total_reads, 2, static_cast<decltype(num_total_reads)>(1000)) << " reads, "
                            << to_si(num_valid_reads, 2, static_cast<decltype(num_total_reads)>(1000)) << " ("
                            << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads) * 100.0
                            << ")% valid reads.";

    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_close(in);
    hts_tpool_destroy(tpool.pool);
}
} // namespace labw::art_modern
