#include "art_modern_config.h"

#include "art_tsam2gsam/lib/cyh_proj_utils.hh"

#include "art_samcmps/lib/AlignmentMiscInfo.hh"
#include "art_samcmps/lib/TranscriptOverlapper.hh"

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

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);

    if (!is_on_mpi_main_process_or_nompi()) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
    init_logger();
    init_file_logger("samgtfcmp");
    print_version("samgtfcmp");
    handle_dumps();
    BOOST_LOG_TRIVIAL(info) << "SYNOPSIS: " << argv[0] << " ref.gtf query.bam out.tsv";
    if (argc != 4) {
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

    samFile* query_sam = CExceptionsProxy::assert_not_null(sam_open(argv[2], get_sam_mode(argv[2], false)),
        USED_HTSLIB_NAME, std::string("Failed to open HTS file ") + argv[2]);
    sam_hdr_t* query_hdr = CExceptionsProxy::assert_not_null(
        sam_hdr_read(query_sam), USED_HTSLIB_NAME, std::string("Failed to get header from HTS file ") + argv[2]);

    auto gtf_stream = std::ifstream(argv[1]);
    TranscriptOverlapper ref_gtf(gtf_stream);
    gtf_stream.close();

    std::ofstream out(argv[3]);
    // clang-format off
    out << "QUERY_QNAME"
        << "\t" << "REF_QNAME"
        << "\t" << "REF_EXONIC_LEN"
        << "\t" << "QUERY_EXONIC_LEN"
        << "\t" << "REF_TOTAL_LEN"
        << "\t" << "QUERY_TOTAL_LEN"
        << "\t" << "REF_N_SS"
        << "\t" << "QUERY_N_SS"
        << "\t" << "REF_N_EXON"
        << "\t" << "QUERY_N_EXON"
        << "\t" << "N_OVERLAPPING_BASE"
        << "\t" << "N_OVERLAPPING_SS"
        << "\t" << "N_OVERLAPPING_EXON"
        << std::endl;
    // clang-format on
    am_readnum_t num_reads = 0;
    bam1_t* query_aln = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");
    int read_retv = 0;
    std::string query_contig_name;

    while (true) {
        read_retv = sam_read1(query_sam, query_hdr, query_aln);
        if (read_retv == SAM_READ_EOF) {
            break;
        }
        if (((query_aln->core.flag & BAM_FSUPPLEMENTARY) | (query_aln->core.flag & BAM_FSECONDARY)) != 0) {
            continue; // TODO: Current version ignores secondary & supplementary alignments.
        }
        query_contig_name = std::string(sam_hdr_tid2name(query_hdr, query_aln->core.tid));

        if ((query_aln->core.flag & BAM_FUNMAP) != 0) {
            continue; // TODO: Current version ignores unmapped alignments.
        }
        AlignmentMiscInfo query_misc;
        query_misc.from_bam_record(query_hdr, query_aln);

        for (const auto* ref_misc :
            ref_gtf.get_overlapping_transcripts(query_contig_name, static_cast<int>(query_aln->core.pos),
                static_cast<int>(query_aln->core.pos
                    + bam_cigar2rlen(static_cast<int>(query_aln->core.n_cigar), bam_get_cigar(query_aln))))) {
            const auto no_base = n_overlapping_base(*ref_misc, query_misc);
            const auto no_ss = n_overlapping_ss(*ref_misc, query_misc);
            const auto no_exons = n_overlapping_exons(*ref_misc, query_misc);
            if (no_base + no_ss + no_exons == 0) {
                continue;
            }
            // clang-format off
        out << query_misc.name
            << "\t"<< ref_misc->name
            << "\t" << ref_misc->exonic_length
            << "\t" << query_misc.exonic_length
            << "\t" << ref_misc->end - ref_misc->start
            << "\t" << bam_cigar2qlen(static_cast<int>(query_aln->core.n_cigar), bam_get_cigar(query_aln))
            << "\t" << ref_misc->ss_starts.size()
            << "\t" << query_misc.ss_starts.size()
            << "\t" << ref_misc->exon_starts.size()
            << "\t" << query_misc.exon_starts.size()
            << "\t" << no_base
            << "\t" << no_ss
            << "\t" << no_exons
            << std::endl;
            // clang-format on
        }

        num_reads++;
        if (num_reads > 1000) {
            std::cerr << "Processed " << num_reads << " reads." << std::endl;
            num_reads = 0;
        }
    }

    sam_hdr_destroy(query_hdr);
    bam_destroy1(query_aln);
    CExceptionsProxy::assert_numeric(sam_close(query_sam), USED_HTSLIB_NAME, "Failed to close SAM/BAM file", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%ws wall, %us user + %ss system = %ts CPU (%p%)");
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
