#include "tsam2gsam/exe/main_fn.hh"

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "tsam2gsam/lib/Transcript.hh"
#include "tsam2gsam/lib/cyh_proj_utils.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <htslib/hts.h>

#include <fmt/format.h>

#include <boost/log/trivial.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#ifdef CEU_CM_IS_DEBUG
#include <ostream>
#endif
#include <string>
#include <vector>

namespace labw::art_modern {
#ifdef CEU_CM_IS_DEBUG

void convert_transcript_to_genome_alignment(
    bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript, std::ostream& cigar_trace)
#else
void convert_transcript_to_genome_alignment(bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript)
#endif
{
    std::string const qname = bam_get_qname(t_aln);
    // Generate OA tag
    auto const oa_tag = fmt::format("{},{},{},{},{},{};", transcript.transcript_id, t_aln->core.pos,
        ((t_aln->core.flag & BAM_FREVERSE) != 0 ? '-' : '+'),
        cigar_arr_to_str(bam_get_cigar(t_aln), t_aln->core.n_cigar), t_aln->core.qual, "");

    if (transcript.is_reverse) {
        if (bam_get_qual(t_aln)[0] != 0xff) {
            reverse(bam_get_qual(t_aln), t_aln->core.l_qseq);
        }
        reverse(bam_get_cigar(t_aln), t_aln->core.n_cigar);

        auto seq_str_tmp = bam_seq_to_str(t_aln);
        revcomp_inplace(seq_str_tmp);
        // Reverse BAM_FLAG
        t_aln->core.flag ^= BAM_FREVERSE;

        bam_set_seq_eqlen(t_aln, seq_str_tmp.c_str());
        t_aln->core.pos = transcript.unspliced_length - bam_endpos(t_aln);
    }

    /** Number of CIGAR at taln **/
    const auto t_n_cigar = t_aln->core.n_cigar;
    /** CIGAR at taln **/
    const am_cigar_t * t_cigar = bam_get_cigar(t_aln);

    auto const qual_str = bam_qual_to_str(t_aln);
    auto const seq_str = bam_seq_to_str(t_aln);

    hts_pos_t g_aln_start = 0;
    auto g_aln_flag = t_aln->core.flag;
    /** Splice site positions oc corresponding transcript **/
    const auto& ssp = transcript.splice_site_positions;

    clear_pe_flag(g_aln_flag);

    std::vector<am_cigar_t> g_cigar;
    g_cigar.reserve(t_n_cigar + ssp.size());

    /** Position on unspliced transcript 5' to 3' **/
    hts_pos_t pos_on_transcript = 0;

    /** Position on read **/
    hts_pos_t pos_on_read = 0;

    /** Position on genome **/
    hts_pos_t pos_on_genome = 0;

    /** CIGAR operation **/
    hts_pos_t this_cigar_ops = 0;

    /** CIGAR length **/
    hts_pos_t this_cigar_len = 0;

    /** CIGAR type. See bam_cigar_type for details. **/
    am_cigar_type_t this_cigar_type = 0;

    // Miscellaneous helper variables
    std::size_t lower_idx = 0;
    std::size_t upper_idx = 0;
    hts_pos_t remaining_exon_length = 0;
    hts_pos_t remaining_cigar_length = 0;

#ifdef CEU_CM_IS_DEBUG
    cigar_trace << ">" << qname << "\n";
#endif

    /** Position on unspliced transcript 5' to 3' **/
    pos_on_transcript = t_aln->core.pos;
    pos_on_read = 0;
    pos_on_genome = transcript.start + t_aln->core.pos;

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
            BOOST_LOG_TRIVIAL(fatal) << "Error: Unknown CIGAR type " << std::to_string(this_cigar_type)
                                     << " for CIGAR operation " << bam_cigar_opchr(this_cigar_ops) << " in read "
                                     << qname << std::endl;
            abort_mpi();
        }
#ifdef CEU_CM_IS_DEBUG
        cigar_trace << "\n";
#endif
    }
    if (pos_on_read != t_aln->core.l_qseq) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: pos_on_read (" << pos_on_read << ") != t_aln->core.l_qseq ("
                                 << t_aln->core.l_qseq << ")" << std::endl;
        abort_mpi();
    }
    g_cigar = merge_cigars(g_cigar);

    // Replace g_cigar
    CExceptionsProxy::assert_numeric(
        bam_set1(
            /** bam **/ g_aln,
            /** l_qname **/ qname.length(),
            /** qname **/ qname.c_str(),
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
            /** l_aux **/ bam_get_l_aux(t_aln) + (/** XT **/ 2 + 1 + transcript.transcript_id.size() + 1)
                + (/** XI **/ 2 + 1 + 4))
            + (/** OA **/ 2 + 1 + oa_tag.size() + 1),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record for " + qname, false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    // Copy existing tags
    std::memcpy(bam_get_aux(g_aln), bam_get_aux(t_aln), bam_get_l_aux(t_aln));
    // Add XT tag
    CExceptionsProxy::assert_numeric(bam_aux_update_str(g_aln, "XT", static_cast<int>(transcript.transcript_id.size()),
                                         transcript.transcript_id.c_str()),
        USED_HTSLIB_NAME, "Failed to add XT tag for " + qname, false, CExceptionsProxy::EXPECTATION::ZERO);
    // Add XI tag
    CExceptionsProxy::assert_numeric(bam_aux_update_int(g_aln, "XI", transcript.tid_on_transcriptome), USED_HTSLIB_NAME,
        "Failed to add XI tag for " + qname, false, CExceptionsProxy::EXPECTATION::ZERO);
    // Add OA tag
    CExceptionsProxy::assert_numeric(bam_aux_update_str(g_aln, "OA", static_cast<int>(oa_tag.size()), oa_tag.c_str()),
        USED_HTSLIB_NAME, "Failed to add OA tag for " + qname, false, CExceptionsProxy::EXPECTATION::ZERO);
}

void populate_ghdr(sam_hdr_t* ghdr, faidx_t* faidx, int argc, char** argv)
{
    CExceptionsProxy::assert_numeric(sam_hdr_remove_lines(ghdr, "SQ", nullptr, nullptr), USED_HTSLIB_NAME,
        "Failed to remove existing SQ lines from SAM header.", false, CExceptionsProxy::EXPECTATION::ZERO);
    for (int i = 0; i < faidx_nseq(faidx); ++i) {
        const auto* name = CExceptionsProxy::assert_not_null(faidx_iseq(faidx, i), USED_HTSLIB_NAME,
            fmt::format("Failed to get sequence name for index {} from FASTA index.", i));
        auto const len = faidx_seq_len64(faidx, name);
        CExceptionsProxy::assert_numeric(
            sam_hdr_add_line(ghdr, "SQ", "SN", name, "LN", std::to_string(len).c_str(), nullptr), USED_HTSLIB_NAME,
            std::string("Failed to populate SQ tag to SAM header for seq ") + name, false,
            CExceptionsProxy::EXPECTATION::ZERO);
    }
    std::string const args = join(std::vector<std::string>(argv, argv + argc), " ");
    CExceptionsProxy::assert_numeric(sam_hdr_change_HD(ghdr, "SO", "unknown"), USED_HTSLIB_NAME,
        "Failed to populate SO tag to SAM header.", false, CExceptionsProxy::EXPECTATION::ZERO);
    CExceptionsProxy::assert_numeric(sam_hdr_add_line(ghdr, "PG", "ID", "tsam2gsam", "PN", "tsam2gsam", "VN",
                                         TSAM2GSAM_VERSION, "CL", args.c_str(), nullptr),
        USED_HTSLIB_NAME, "Failed to populate PG tag to SAM header.", false, CExceptionsProxy::EXPECTATION::ZERO);
}
} // namespace labw::art_modern
