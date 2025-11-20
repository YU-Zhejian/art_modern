/**
 * This program tries to deal with the situation where
 * both aligners do not report sequences that are failed to align.
 *
 * Comparing 2 alignment files on real data with 1 unaligned BAM generated from FASTQ files.
 *
 */
#include "art_modern_config.h"

#include "art_tsam2gsam/lib/cyh_proj_utils.hh"

#include "art_samcmps/lib/AlignmentMiscInfo.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/dump_utils.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

// Boost timer
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <htslib/sam.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace labw::art_modern;

namespace {
inline bool is_valid(const bam1_t* aln)
{
    return ((((aln->core.flag & BAM_FSUPPLEMENTARY) | (aln->core.flag & BAM_FSECONDARY) | (aln->core.flag & BAM_FUNMAP))
        != 0));
}
} // namespace

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);

    if (!is_on_mpi_main_process_or_nompi()) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
    init_logger();
    init_file_logger("samcmp3");
    print_version("samcmp3");
    handle_dumps();
    BOOST_LOG_TRIVIAL(info) << "SYNOPSIS: " << argv[0] << " ref.bam query1.bam query2.bam out.tsv";
    if (argc != 5) {
        BOOST_LOG_TRIVIAL(fatal) << "ERROR: Invalid number of arguments!";
        abort_mpi();
    }

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    samFile* ref_sam = CExceptionsProxy::assert_not_null(sam_open(argv[1], get_sam_mode(argv[1], false)),
        USED_HTSLIB_NAME, std::string("Failed to open HTS file ") + argv[1]);
    samFile* query1_sam = CExceptionsProxy::assert_not_null(sam_open(argv[2], get_sam_mode(argv[2], false)),
        USED_HTSLIB_NAME, std::string("Failed to open HTS file ") + argv[2]);
    samFile* query2_sam = CExceptionsProxy::assert_not_null(sam_open(argv[3], get_sam_mode(argv[3], false)),
        USED_HTSLIB_NAME, std::string("Failed to open HTS file ") + argv[3]);
    std::ofstream out(argv[4]);
    // clang-format off
    out << "REF_QNAME"
        << "\t" << "QUERY1_EXONIC_LEN"
        << "\t" << "QUERY2_EXONIC_LEN"
        << "\t" << "QUERY1_TOTAL_LEN"
        << "\t" << "QUERY2_TOTAL_LEN"
        << "\t" << "QUERY1_N_SS"
        << "\t" << "QUERY2_N_SS"
        << "\t" << "QUERY1_N_EXON"
        << "\t" << "QUERY2_N_EXON"
        << "\t" << "N_OVERLAPPING_BASE"
        << "\t" << "N_OVERLAPPING_SS"
        << "\t" << "N_OVERLAPPING_EXON"
        << std::endl;
    // clang-format on
    am_readnum_t num_reads = 0;
    am_readnum_t total_num_reads = 0;
    sam_hdr_t* ref_hdr = CExceptionsProxy::assert_not_null(
        sam_hdr_read(ref_sam), USED_HTSLIB_NAME, std::string("Failed to get header from HTS file ") + argv[1]);
    sam_hdr_t* query1_hdr = CExceptionsProxy::assert_not_null(
        sam_hdr_read(query1_sam), USED_HTSLIB_NAME, std::string("Failed to get header from HTS file ") + argv[2]);
    sam_hdr_t* query2_hdr = CExceptionsProxy::assert_not_null(
        sam_hdr_read(query2_sam), USED_HTSLIB_NAME, std::string("Failed to get header from HTS file ") + argv[3]);
    bam1_t* ref_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    bam1_t* query1_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    bam1_t* query2_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    AlignmentMiscInfo query1_misc;
    AlignmentMiscInfo query2_misc;
    int read_retv = 0;
    bool eof_encountered = false;
    std::string ref_qname;
    std::string query1_qname;
    std::string query2_qname;
    // Load the 1st query from each library
    do {
        read_retv = sam_read1(query1_sam, query1_hdr, query1_aln);
        if (read_retv == SAM_READ_EOF) {
            eof_encountered = true;
            break;
        }
    } while (is_valid(query1_aln));

    do {
        read_retv = sam_read1(query2_sam, query2_hdr, query2_aln);
        if (read_retv == SAM_READ_EOF) {
            eof_encountered = true;
            break;
        }
    } while (is_valid(query2_aln));

    while (!eof_encountered) {
        query1_qname = bam_get_qname(query1_aln);
        query2_qname = bam_get_qname(query2_aln);

        read_retv = sam_read1(ref_sam, ref_hdr, ref_aln);
        if (read_retv == SAM_READ_EOF) {
            eof_encountered = true;
            break;
        }
        ref_qname = bam_get_qname(ref_aln);
        if (ref_qname == query1_qname) {
            query1_misc.from_bam_record(query1_hdr, query1_aln);
        } else {
            query1_misc.reset();
        }
        if (ref_qname == query2_qname) {
            query2_misc.from_bam_record(query2_hdr, query2_aln);
        } else {
            query2_misc.reset();
        }

        // clang-format off
        out << ref_qname
            << "\t" << query1_misc.exonic_length // Default to 0
            << "\t" << query2_misc.exonic_length // Default to 0
            << "\t" << (query1_misc.is_unaligned()? 0 :bam_cigar2qlen(static_cast<int>(query1_aln->core.n_cigar), bam_get_cigar(query1_aln)))
            << "\t" << (query2_misc.is_unaligned()? 0 :bam_cigar2qlen(static_cast<int>(query2_aln->core.n_cigar), bam_get_cigar(query2_aln)))
            << "\t" << query1_misc.ss_starts.size() // Default to 0
            << "\t" << query2_misc.ss_starts.size() // Default to 0
            << "\t" << query1_misc.exon_starts.size() // Default to 0
            << "\t" << query2_misc.exon_starts.size() // Default to 0
            << "\t" << n_overlapping_base(query1_misc, query2_misc)
            << "\t" << n_overlapping_ss(query1_misc, query2_misc)
            << "\t" << n_overlapping_exons(query1_misc, query2_misc)
            << std::endl;
        // clang-format on

        if (ref_qname == query1_qname) {
            do {
                read_retv = sam_read1(query1_sam, query1_hdr, query1_aln);
                if (read_retv == SAM_READ_EOF) {
                    eof_encountered = true;
                    break;
                }
            } while (is_valid(query1_aln));
        }
        if (ref_qname == query2_qname) {
            do {
                read_retv = sam_read1(query2_sam, query2_hdr, query2_aln);
                if (read_retv == SAM_READ_EOF) {
                    eof_encountered = true;
                    break;
                }
            } while (is_valid(query2_aln));
        }

        num_reads++;
        total_num_reads++;
        if (num_reads > 1000) {
            std::cerr << "Processed " << num_reads << " reads." << std::endl;
            num_reads = 0;
        }
    }
    std::cerr << "# total reads: " << total_num_reads << std::endl;

    sam_hdr_destroy(ref_hdr);
    sam_hdr_destroy(query1_hdr);
    sam_hdr_destroy(query2_hdr);
    bam_destroy1(ref_aln);
    bam_destroy1(query1_aln);
    bam_destroy1(query2_aln);
    CExceptionsProxy::assert_numeric(sam_close(ref_sam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(sam_close(query1_sam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(sam_close(query2_sam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%ws wall, %us user + %ss system = %ts CPU (%p%)");
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
