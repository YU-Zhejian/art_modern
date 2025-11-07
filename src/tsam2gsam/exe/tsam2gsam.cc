/**
 * @file tsam2gsam.cc
 * @author YU Zhejian (zhejianyu@intl.zju.edu.cn)
 * @brief Convert transcriptome-aligned SAM to genome-aligned SAM
 * @version 0.1
 * @date 2024-11-09
 *
 * @copyright Copyright (c) 2024 YU Zhejian
 *
 */

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "tsam2gsam/lib/Transcript.hh"
#include "tsam2gsam/lib/cyh_proj_utils.hh"
#include "tsam2gsam/lib/gffread_bed_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/utils/seq_utils.hh"
#include "libam_support/utils/dump_utils.hh"

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace labw::art_modern;
namespace {

/**
 * @brief Convert transcript-aligned read to genome-aligned read.
 *
 * @code
 *
 *          Exon 0         Exon 1                    Exon 2
 *
 *  TALN   |<====|--------|====|--------------------|=============|
 *  GALN       |=|--------|====|--------------------|==>|
 *
 * @endcode
 *
 * @param t_aln Input transcript-aligned read.
 * @param g_aln Output genome-aligned read.
 * @param transcript The transcript.
 * @return true if there's no error; false otherwise
 */
#ifdef CEU_CM_IS_DEBUG
void convert_transcript_to_genome_alignment(
    bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript, std::ostream& cigar_trace)
#else
void convert_transcript_to_genome_alignment(bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript)
#endif
{
    // Generate OA tag
    std::ostringstream oa_ss;
    oa_ss << /** Contig **/ transcript.transcript_id << ',' << /** Pos **/ t_aln->core.pos << ','
          << /** Strand **/ ((t_aln->core.flag & BAM_FREVERSE) != 0 ? '-' : '+') << ','
          << /** CIGAR **/ cigar_arr_to_str(bam_get_cigar(t_aln), t_aln->core.n_cigar) << ','
          << /** MAPQ **/ std::to_string(t_aln->core.qual) << ',' << /** NM **/ "" << ';';
    auto const oa_tag = oa_ss.str();

    if (transcript.is_reverse) {
        reverse(bam_get_qual(t_aln), t_aln->core.l_qseq);
        reverse(bam_get_cigar(t_aln), t_aln->core.n_cigar);

        auto seq_str = bam_seq_to_str(t_aln);
        revcomp_inplace(seq_str);
        // Reverse BAM_FLAG
        t_aln->core.flag ^= BAM_FREVERSE;

        bam_set_seq_eqlen(t_aln, seq_str.c_str());
        t_aln->core.pos = transcript.unspliced_length - bam_endpos(t_aln);
    }

    /** Number of CIGAR at taln **/
    const uint32_t t_n_cigar = t_aln->core.n_cigar;
    /** CIGAR at taln **/
    const uint32_t* t_cigar = bam_get_cigar(t_aln);

    auto qual_str = bam_qual_to_str(t_aln);
    auto seq_str = bam_seq_to_str(t_aln);

    int32_t g_aln_start = 0;
    auto g_aln_flag = t_aln->core.flag;
    /** Splice site positions oc corresponding transcript **/
    const auto& ssp = transcript.splice_site_positions;

    clear_pe_flag(g_aln_flag);

    std::vector<uint32_t> g_cigar;
    g_cigar.reserve(t_n_cigar + ssp.size());

    /** Position on unspliced transcript 5' to 3' **/
    int32_t pos_on_transcript = 0;

    /** Position on read **/
    int32_t pos_on_read = 0;

    /** Position on genome **/
    int32_t pos_on_genome = 0;

    /** CIGAR operation **/
    int32_t this_cigar_ops = 0;

    /** CIGAR length **/
    int32_t this_cigar_len = 0;

    /** CIGAR type. See bam_cigar_type for details. **/
    int this_cigar_type = 0;

    // Miscellaneous helper variables
    std::size_t lower_idx = 0;
    std::size_t upper_idx = 0;
    int32_t remaining_exon_length = 0;
    int32_t remaining_cigar_length = 0;

#ifdef CEU_CM_IS_DEBUG
    cigar_trace << ">" << bam_get_qname(t_aln) << "\n";
#endif

    /** Position on unspliced transcript 5' to 3' **/
    pos_on_transcript = static_cast<int32_t>(t_aln->core.pos);
    pos_on_read = 0;
    pos_on_genome = static_cast<int32_t>(transcript.start + t_aln->core.pos);

    // Get alignment start by populating introns
    lower_idx = std::lower_bound(ssp.begin(), ssp.end(), /** Transcript start **/ 0) - ssp.begin();
    upper_idx = std::upper_bound(ssp.begin(), ssp.end(), pos_on_transcript - 1) - ssp.begin();
    for (auto ss_idx = lower_idx; ss_idx < upper_idx; ++ss_idx) {
        pos_on_genome += transcript.splice_site_lengths[ss_idx];
    }
    g_aln_start = pos_on_genome;

    for (decltype(t_aln->core.n_cigar) t_cigar_idx = 0; t_cigar_idx < t_n_cigar; ++t_cigar_idx) {
        this_cigar_ops = bam_cigar_op(t_cigar[t_cigar_idx]);
        this_cigar_len = bam_cigar_oplen(t_cigar[t_cigar_idx]);
        this_cigar_type = bam_cigar_type(this_cigar_ops);
#ifdef CEU_CM_IS_DEBUG
        cigar_trace << this_cigar_len << bam_cigar_opchr(this_cigar_ops) << " ";
#endif
        switch (this_cigar_type) {
        case CONSUME_QUERY_AND_REFERENCE:
            // If this CIGAR segment contains an intron, make it 2 cigar segments with N inserted.
            // Note that it is possible to introduce more than 1 intron.
            lower_idx = std::lower_bound(ssp.begin(), ssp.end(), pos_on_transcript) - ssp.begin();
            upper_idx = std::upper_bound(ssp.begin(), ssp.end(), pos_on_transcript + this_cigar_len - 1) - ssp.begin();
            remaining_cigar_length = this_cigar_len;
            for (auto ss_idx = lower_idx; ss_idx < upper_idx; ++ss_idx) {
                remaining_exon_length = ssp[ss_idx] - pos_on_transcript;
                remaining_cigar_length -= remaining_exon_length;
                g_cigar.emplace_back(bam_cigar_gen(remaining_exon_length, this_cigar_ops));
#ifdef CEU_CM_IS_DEBUG
                cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
                // Skip exonic regions on the genome
                pos_on_transcript += remaining_exon_length;
                pos_on_genome += remaining_exon_length;
                pos_on_read += remaining_exon_length;

                // Skip intronic regions on the genome
                pos_on_genome += transcript.splice_site_lengths[ss_idx];
                g_cigar.emplace_back(bam_cigar_gen(transcript.splice_site_lengths[ss_idx], BAM_CREF_SKIP));
#ifdef CEU_CM_IS_DEBUG
                cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
            }
            if (remaining_cigar_length > 0) {
                remaining_exon_length = remaining_cigar_length;
                g_cigar.emplace_back(bam_cigar_gen(remaining_exon_length, this_cigar_ops));
#ifdef CEU_CM_IS_DEBUG
                cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
                pos_on_transcript += remaining_exon_length;
                pos_on_genome += remaining_exon_length;
                pos_on_read += remaining_exon_length;
            }
            break;
        case CONSUME_NEITHER_QUERY_NOR_REFERENCE:
            g_cigar.emplace_back(bam_cigar_gen(this_cigar_len, this_cigar_ops));
#ifdef CEU_CM_IS_DEBUG
            cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
            // Do nothing!
            break;
        case CONSUME_QUERY:
            pos_on_read += this_cigar_len;
            g_cigar.emplace_back(bam_cigar_gen(this_cigar_len, this_cigar_ops));
#ifdef CEU_CM_IS_DEBUG
            cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
            // No advancement on pos_on_transcript and pos_on_genome
            break;
        case CONSUME_REFERENCE:
            // If the deletion & refskip overlaps with an intron, skip it!
            // Note that it is possible to introduce more than 1 intron.
            // Get splice sites within the cigar range.
            lower_idx = std::lower_bound(ssp.begin(), ssp.end(), pos_on_transcript) - ssp.begin();
            upper_idx = std::upper_bound(ssp.begin(), ssp.end(), pos_on_transcript + this_cigar_len - 1) - ssp.begin();

            pos_on_transcript += this_cigar_len; // Skip exonic regions on the transcript
            pos_on_genome += this_cigar_len; // Skip exonic regions on the genome

            // Skip intronic regions on the genome
            for (auto ss_idx = lower_idx; ss_idx < upper_idx; ++ss_idx) {
                this_cigar_ops = BAM_CREF_SKIP; // Deletion becomes N
                pos_on_genome += transcript.splice_site_lengths[ss_idx];
                this_cigar_len += transcript.splice_site_lengths[ss_idx];
            }
            g_cigar.emplace_back(bam_cigar_gen(this_cigar_len, this_cigar_ops));
#ifdef CEU_CM_IS_DEBUG
            cigar_trace << cigar_arr_to_str(&g_cigar.back(), 1);
#endif
            // No advancement on pos_of_read
            break;
        default: // Error
            std::abort();
        }
#ifdef CEU_CM_IS_DEBUG
        cigar_trace << "\n";
#endif
    }
    if (pos_on_read != t_aln->core.l_qseq) {
        std::cerr << "Error: pos_on_read (" << pos_on_read << ") != t_aln->core.l_qseq (" << t_aln->core.l_qseq << ")"
                  << std::endl;
        std::abort();
    }
    g_cigar = merge_cigars(g_cigar);

    // Replace g_cigar
    require_not_negative(bam_set1(
                             /** bam **/ g_aln,
                             /** l_qname **/ std::strlen(bam_get_qname(t_aln)),
                             /** qname **/ bam_get_qname(t_aln),
                             /** flag **/ g_aln_flag,
                             /** tid **/ transcript.tid_on_genome,
                             /** pos **/ g_aln_start,
                             /** mapq **/ t_aln->core.qual, // Ridiculous
                             /** n_cigar **/ static_cast<uint32_t>(g_cigar.size()),
                             /** cigar **/ g_cigar.data(),
                             /** mtid **/ 0,
                             /** mpos **/ 0,
                             /** isize **/ 0,
                             /** l_seq **/ t_aln->core.l_qseq,
                             /** seq **/ seq_str.c_str(),
                             /** qual **/ qual_str.c_str(),
                             /** l_aux **/ bam_get_l_aux(t_aln)
                                 + (/** XT **/ 2 + 1 + transcript.transcript_id.size() + 1) + (/** XI **/ 2 + 1 + 4))
        + (/** OA **/ 2 + 1 + oa_tag.size() + 1));

    // Copy existing tags
    std::memcpy(bam_get_aux(g_aln), bam_get_aux(t_aln), bam_get_l_aux(t_aln));
    // Add XT tag
    require_zero(bam_aux_update_str(
        g_aln, "XT", static_cast<int>(transcript.transcript_id.size()), transcript.transcript_id.c_str()));
    // Add XI tag
    require_zero(bam_aux_update_int(g_aln, "XI", transcript.tid_on_transcriptome));
    // Add OA tag
    require_zero(bam_aux_update_str(g_aln, "OA", static_cast<int>(oa_tag.size()), oa_tag.c_str()));
}

void populate_ghdr(sam_hdr_t* ghdr, faidx_t* faidx, int argc, char const* argv[])
{
    require_zero(sam_hdr_remove_lines(ghdr, "SQ", nullptr, nullptr));
    for (int i = 0; i < faidx_nseq(faidx); ++i) {
        const auto* name = require_not_null(faidx_iseq(faidx, i));
        int const len = faidx_seq_len(faidx, name);
        require_zero(sam_hdr_add_line(ghdr, "SQ", "SN", name, "LN", std::to_string(len).c_str(), nullptr));
    }
    std::stringstream args_ss;
    for (int i = 0; i < argc - 1; ++i) {
        args_ss << argv[i] << " ";
    }
    args_ss << argv[argc - 1];
    require_zero(sam_hdr_change_HD(ghdr, "SO", "unknown"));
    require_zero(sam_hdr_add_line(
            ghdr, "PG", "ID", "tsam2gsam", "PN", "tsam2gsam", "VN", TSAM2GSAM_VERSION, "CL", args_ss.str().c_str(), nullptr));
}

}; // namespace
int main(int argc, char const* argv[])
{
    print_version("tsam2gsam");
    handle_dumps();
    std::cout << "SYNOPSIS: " << argv[0]
              << " taln.bam galn.bam reference_annotation.gtf.gffread.bed reference_assembly.fa" << std::endl;
    if (argc != 5) {
        std::cerr << "ERROR: Invalid arguments." << std::endl;
        return EXIT_FAILURE;
    }
    samFile* tsam = require_not_null(sam_open(argv[1], get_sam_mode(argv[1], false)));
    sam_hdr_t* thdr = require_not_null(sam_hdr_read(tsam));
    samFile* gsam = require_not_null(sam_open(argv[2], get_sam_mode(argv[1], true)));
    hts_set_threads(gsam, N_THREADS_FOR_HTSLIB); // 10 concurrent threads for compression
#ifdef CEU_CM_IS_DEBUG
    std::ofstream cigar_trace(std::string(argv[2]) + ".ctrace");
#endif

    // Formulate SAM header for GSAM
    std::cerr << "INFO: Formulating SAM header..." << std::endl;
    faidx_t* faidx = require_not_null(fai_load(argv[4]));
    sam_hdr_t* ghdr = require_not_null(sam_hdr_dup(thdr));
    populate_ghdr(ghdr, faidx, argc, argv);
    require_zero(sam_hdr_write(gsam, ghdr));
    fai_destroy(faidx);
    std::cerr << "INFO: SAM header written." << std::endl;

    // Read GTF
    std::cerr << "INFO: Reading GTF..." << std::endl;
    std::ifstream gtf_stream(argv[3]);
    auto transcript_id_to_transcript_map = read_gffutils_bed(thdr, ghdr, gtf_stream);
    gtf_stream.close();
    std::cerr << "INFO: " << transcript_id_to_transcript_map.size() << " transcripts read." << std::endl;

    // Start the main loop
    std::cerr << "INFO: Start main loop..." << std::endl;
    bam1_t* t_aln = require_not_null(bam_init1());
    bam1_t* g_aln = require_not_null(bam_init1());
    std::size_t num_correct_transcripts = 0;

    while (sam_read1(tsam, thdr, t_aln) >= 0) {
        // Passthrough of unmapped transcripts
        if ((t_aln->core.flag & BAM_FUNMAP) != 0) {
            num_correct_transcripts++;
            require_not_negative(sam_write1(gsam, ghdr, t_aln)); // Write unmapped reads as-is.
        }
        auto transcript_name = std::string(require_not_null(sam_hdr_tid2name(thdr, t_aln->core.tid)));
        if (transcript_id_to_transcript_map.find(transcript_name)
            == transcript_id_to_transcript_map.end()) {
            std::cerr << "ERROR: Transcript " << transcript_name << " not found in GTF." << std::endl;
            return EXIT_FAILURE;
        }
        convert_transcript_to_genome_alignment(t_aln, g_aln, transcript_id_to_transcript_map.at(transcript_name)
#ifdef CEU_CM_IS_DEBUG
                                                                 ,
            cigar_trace
#endif
        );
        require_not_negative(sam_write1(gsam, ghdr, g_aln));
        num_correct_transcripts++;
        if (num_correct_transcripts % 10000 == 0) {
            std::cerr << "INFO: " << num_correct_transcripts << " transcripts processed." << std::endl;
        }
    }
    bam_destroy1(t_aln);
    bam_destroy1(g_aln);
    std::cerr << "INFO: " << num_correct_transcripts << " transcripts were mapped correctly." << std::endl;

    // Cleanup
    std::cerr << "INFO: Cleaning up..." << std::endl;
    sam_hdr_destroy(thdr);
    sam_hdr_destroy(ghdr);
    require_zero(sam_close(tsam));
    require_zero(sam_close(gsam));
#ifdef CEU_CM_IS_DEBUG
    cigar_trace.close();
#endif
    std::cerr << "INFO: Done." << std::endl;
    return EXIT_SUCCESS;
}
