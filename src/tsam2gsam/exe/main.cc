/**
 * @file tsam2gsam.cc
 * @author YU Zhejian (zhejianyu@intl.zju.edu.cn)
 * @brief Convert transcriptome-aligned SAM to genome-aligned SAM
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024 YU Zhejian
 *
 * TODO: This code sufferes from QUAL string issues.
 */

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "tsam2gsam/exe/main_fn.hh"

#include "tsam2gsam/lib/cyh_proj_utils.hh"
#include "tsam2gsam/lib/gffread_bed_utils.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/utils/dump_utils.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

// Boost timer
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

using namespace labw::art_modern;

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);

    if (!is_on_mpi_main_process_or_nompi()) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
    init_logger();
    init_file_logger("tsam2gsam");
    print_version("tsam2gsam");
    handle_dumps();
    BOOST_LOG_TRIVIAL(info) << "SYNOPSIS: " << argv[0]
                            << " taln.bam galn.bam reference_annotation.gtf.gffread.bed reference_assembly.fa";
    if (argc != 5) {
        BOOST_LOG_TRIVIAL(fatal) << "Invalid arguments.";
        abort_mpi();
    }

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif

    samFile* tsam = CExceptionsProxy::assert_not_null(sam_open(argv[1], get_sam_mode(argv[1], false)), USED_HTSLIB_NAME,
        std::string("Failed to open HTS file ") + argv[1]);
    sam_hdr_t* thdr = CExceptionsProxy::assert_not_null(
        sam_hdr_read(tsam), USED_HTSLIB_NAME, std::string("Failed to get header from HTS file ") + argv[1]);
    samFile* gsam = CExceptionsProxy::assert_not_null(sam_open(argv[2], get_sam_mode(argv[1], true)), USED_HTSLIB_NAME,
        std::string("Failed to open HTS file ") + argv[2]);
    hts_set_threads(gsam, N_THREADS_FOR_HTSLIB); // 10 concurrent threads for compression
#ifdef CEU_CM_IS_DEBUG
    std::ofstream cigar_trace(std::string(argv[2]) + ".ctrace");
#endif

    // Formulate SAM header for GSAM
    BOOST_LOG_TRIVIAL(info) << "Formulating SAM header...";
    faidx_t* faidx
        = CExceptionsProxy::assert_not_null(fai_load(argv[4]), USED_HTSLIB_NAME, "Failed to load FASTA index.");
    sam_hdr_t* ghdr
        = CExceptionsProxy::assert_not_null(sam_hdr_dup(thdr), USED_HTSLIB_NAME, "Failed to duplicate SAM header.");
    populate_ghdr(ghdr, faidx, argc, argv);
    CExceptionsProxy::assert_numeric(sam_hdr_write(gsam, ghdr), USED_HTSLIB_NAME,
        "Failed to write SAM header to output file.", false, CExceptionsProxy::EXPECTATION::ZERO);
    fai_destroy(faidx);
    BOOST_LOG_TRIVIAL(info) << "SAM header written.";

    // Read GTF
    BOOST_LOG_TRIVIAL(info) << "Reading GTF...";
    std::ifstream gtf_stream(argv[3]);
    const auto& transcript_id_to_transcript_map = read_gffutils_bed(thdr, ghdr, gtf_stream);
    gtf_stream.close();
    BOOST_LOG_TRIVIAL(info) << transcript_id_to_transcript_map.size() << " transcripts read.";

    // Start the main loop
    BOOST_LOG_TRIVIAL(info) << "Start main loop...";
    bam1_t* t_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    bam1_t* g_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    std::size_t num_correct_transcripts = 0;

    while (sam_read1(tsam, thdr, t_aln) >= 0) {
        // Passthrough of unmapped transcripts
        if ((t_aln->core.flag & BAM_FUNMAP) != 0) {
            num_correct_transcripts++;
            // Write unmapped reads as-is.
            CExceptionsProxy::assert_numeric(sam_write1(gsam, ghdr, t_aln), USED_HTSLIB_NAME,
                std::string("Failed to write SAM entry '") + bam_get_qname(t_aln) + "' to output file.", false,
                CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        }
        const auto transcript_name = std::string(CExceptionsProxy::assert_not_null(
            sam_hdr_tid2name(thdr, t_aln->core.tid), USED_HTSLIB_NAME, "Failed to get transcript name."));
        const auto transcript_iter = transcript_id_to_transcript_map.find(transcript_name);
        if (transcript_iter == transcript_id_to_transcript_map.end()) {
            BOOST_LOG_TRIVIAL(error) << "Transcript " << transcript_name << " not found in GTF, skipping...";
            continue;
        }
        convert_transcript_to_genome_alignment
#ifdef CEU_CM_IS_DEBUG
            (t_aln, g_aln, transcript_iter->second, cigar_trace);
#else
            (t_aln, g_aln, transcript_iter->second);
#endif
        CExceptionsProxy::assert_numeric(sam_write1(gsam, ghdr, t_aln), USED_HTSLIB_NAME,
            std::string("Failed to write SAM entry '") + bam_get_qname(t_aln) + "' to output file.", false,
            CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        num_correct_transcripts++;
        if (num_correct_transcripts % 10000 == 0) {
            BOOST_LOG_TRIVIAL(info) << num_correct_transcripts << " transcripts processed.";
        }
    }
    bam_destroy1(t_aln);
    bam_destroy1(g_aln);
    BOOST_LOG_TRIVIAL(info) << num_correct_transcripts << " transcripts were mapped correctly.";

    // Cleanup
    BOOST_LOG_TRIVIAL(info) << "Cleaning up...";
    sam_hdr_destroy(thdr);
    sam_hdr_destroy(ghdr);
    CExceptionsProxy::assert_numeric(sam_close(tsam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(sam_close(gsam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
#ifdef CEU_CM_IS_DEBUG
    cigar_trace.close();
#endif

#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%ws wall, %us user + %ss system = %ts CPU (%p%)");
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
