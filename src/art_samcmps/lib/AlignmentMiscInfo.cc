#include "art_samcmps/lib/AlignmentMiscInfo.hh"

#include "art_tsam2gsam/lib/cyh_proj_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/arithmetic_utils.hh"

#include <IITree.h>

#include <ewah/ewah.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <array>
#include <cstdlib>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
void AlignmentMiscInfo::from_bam_record(const bam_hdr_t* hdr, const bam1_t* aln)
{
    reset();
    name = bam_get_qname(aln);
    am_cigar_ops_t this_cigar_ops = 0;
    am_cigar_len_t this_cigar_len = 0;
    am_cigar_type_t this_cigar_type = 0;

    auto pos_on_genome = aln->core.pos;
    start = pos_on_genome;
    [[maybe_unused]] hts_pos_t pos_on_read = 0;
    contig_name = sam_hdr_tid2name(hdr, aln->core.tid);

    const std::size_t n_cigar = aln->core.n_cigar;
    const am_cigar_t* cigar = bam_get_cigar(aln);
    hts_pos_t current_exon_start = -1;
    for (std::size_t cigar_idx = 0; cigar_idx < n_cigar; ++cigar_idx) {
        this_cigar_ops = bam_cigar_op(cigar[cigar_idx]);
        this_cigar_len = static_cast<am_cigar_len_t>(bam_cigar_oplen(cigar[cigar_idx]));
        this_cigar_type = bam_cigar_type(this_cigar_ops);

        switch (this_cigar_type) {
        case CONSUME_QUERY_AND_REFERENCE:
            if (current_exon_start == -1) {
                current_exon_start = pos_on_genome;
            }
            exonic_length += this_cigar_len;
            for (decltype(pos_on_genome) i = pos_on_genome; i < pos_on_genome + this_cigar_len; ++i) {
                base_set.set(i);
            }
            pos_on_genome += this_cigar_len;
            pos_on_read += this_cigar_len;
            break;
        case CONSUME_NEITHER_QUERY_NOR_REFERENCE:
            break;
        case CONSUME_QUERY:
            if (this_cigar_ops == BAM_CINS) {
                // FIXME: emplace_back ERROR! Why?
                insertion_pos.push_back(pos_on_genome);
                insertion_lengths.push_back(this_cigar_len);
            }
            pos_on_read += this_cigar_len;
            break;
        case CONSUME_REFERENCE:
            if (this_cigar_ops == BAM_CDEL) {
                deletion_starts.push_back(pos_on_genome);
                deletion_ends.push_back(pos_on_genome + this_cigar_len);
            } else if (this_cigar_ops == BAM_CREF_SKIP) {
                ss_starts.emplace_back(pos_on_genome);
                ss_ends.emplace_back(pos_on_genome + this_cigar_len);
                exon_starts.emplace_back(current_exon_start);
                exon_ends.emplace_back(pos_on_genome);
                current_exon_start = -1;
            }
            pos_on_genome += this_cigar_len;
            break;
        default: // Error!
            std::abort();
        }
    }
    if (current_exon_start != -1) {
        exon_starts.emplace_back(current_exon_start);
        exon_ends.emplace_back(pos_on_genome);
    }
#ifdef CEU_CM_IS_DEBUG
    if (pos_on_read != bam_cigar2qlen(static_cast<int>(n_cigar), cigar)) {
        std::cerr << "ERROR: Cigar string does not match query length!" << std::endl;
        std::abort();
    }
    if (pos_on_genome != bam_endpos(aln)) {
        std::cerr << "ERROR: Cigar string does not match reference length!" << std::endl;
        std::abort();
    }
#endif
    end = pos_on_genome;
    update_cr_index();
}

void AlignmentMiscInfo::from_gffutils_bed_line(const std::string& line)
{
    reset();
    std::array<std::string, NUM_GFFUTILS_BED_FIELDS> tokens;
    std::string token {};
    int n_exons = 0;
    std::istringstream line_ss(line);
    for (int i = 0; i < NUM_GFFUTILS_BED_FIELDS; ++i) {
        std::getline(line_ss, token, '\t');
        tokens[i] = std::move(token);
    }
    contig_name = std::move(tokens[0]);
    start = std::stoi(tokens[1]);
    end = std::stoi(tokens[2]);
    name = std::move(tokens[3]);

    n_exons = std::stoi(tokens[9]);
    std::istringstream elen_ss(tokens[10]);
    std::istringstream epos_ss(tokens[11]);
    hts_pos_t exon_start = 0;
    hts_pos_t exon_end = 0;
    for (decltype(n_exons) exon_id = 0; exon_id < n_exons; ++exon_id) {
        std::getline(epos_ss, token, ',');
        exon_start = std::stoi(token) + start;
        exon_starts.emplace_back(exon_start);

        std::getline(elen_ss, token, ',');
        exon_end = start + std::stoi(token);
        exon_ends.emplace_back(exon_end);
        exonic_length += exon_end - exon_start;
        for (auto i = exon_start; i < exon_end; ++i) {
            base_set.set(i);
        }
        if (exon_id != 0) {
            ss_starts.emplace_back(exon_ends[exon_id - 1]);
            ss_ends.emplace_back(exon_start);
        }
    }
#ifdef CEU_CM_IS_DEBUG
    if (exon_ends.back() != end) {
        std::abort();
    }
#endif
    update_cr_index();
}
void AlignmentMiscInfo::update_cr_index()
{
    for (decltype(ss_starts.size()) i = 0; i < ss_starts.size(); ++i) {
        ss.add(ss_starts[i], ss_ends[i], i);
        ss_with_err.add(
            am_max(static_cast<hts_pos_t>(0), ss_starts[i] - ALLOWED_ERROR_BASES), ss_ends[i] + ALLOWED_ERROR_BASES, i);
    }
    for (decltype(deletion_starts.size()) i = 0; i < deletion_starts.size(); ++i) {
        deletions.add(deletion_starts[i], deletion_ends[i], i);
        deletions_with_err.add(am_max(static_cast<hts_pos_t>(0), deletion_starts[i] - ALLOWED_ERROR_BASES),
            deletion_ends[i] + ALLOWED_ERROR_BASES, i);
    }
    for (decltype(exon_starts.size()) i = 0; i < exon_starts.size(); ++i) {
        exons.add(exon_starts[i], exon_ends[i], i);
        exons_with_err.add(am_max(static_cast<hts_pos_t>(0), exon_starts[i] - ALLOWED_ERROR_BASES),
            exon_ends[i] + ALLOWED_ERROR_BASES, i);
    }
    ss.index();
    deletions.index();
    ss_with_err.index();
    deletions_with_err.index();
    exons.index();
    exons_with_err.index();
}
bool AlignmentMiscInfo::is_unaligned() const { return contig_name.empty(); }
void AlignmentMiscInfo::reset()
{
    contig_name.clear();
    start = 0;
    end = 0;
    ss_starts.clear();
    ss_ends.clear();
    deletion_starts.clear();
    deletion_ends.clear();
    insertion_pos.clear();
    insertion_lengths.clear();
    exon_starts.clear();
    exon_ends.clear();
    exonic_length = 0;
    base_set = ewah::EWAHBoolArray<std::size_t>();
    ss = IITree<hts_pos_t, std::size_t>();
    deletions = IITree<hts_pos_t, std::size_t>();
    exons = IITree<hts_pos_t, std::size_t>();
    ss_with_err = IITree<hts_pos_t, std::size_t>();
    deletions_with_err = IITree<hts_pos_t, std::size_t>();
    exons_with_err = IITree<hts_pos_t, std::size_t>();
}

std::size_t n_overlapping_base(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query)
{
    if (ref.is_unaligned() || query.is_unaligned() || ref.contig_name != query.contig_name) {
        return 0;
    }
    auto base_set = ref.base_set & query.base_set;
    return base_set.numberOfOnes();
}
std::size_t n_overlapping_ss(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query)
{
    if (ref.is_unaligned() || query.is_unaligned() || ref.contig_name != query.contig_name || ref.ss_starts.empty()
        || query.ss_starts.empty()) {
        return 0;
    }
    std::vector<std::size_t> hit_ids {};

    std::set<std::size_t> ss_matches;

    for (decltype(query.ss_starts.size()) ss_idx_on_query = 0; ss_idx_on_query < query.ss_starts.size();
        ++ss_idx_on_query) {
        ref.ss_with_err.overlap(query.ss_starts[ss_idx_on_query], query.ss_ends[ss_idx_on_query], hit_ids);
        for (const auto& ss_idx_on_ref : hit_ids) {
            if (std::abs(ref.ss_starts[ss_idx_on_ref] - query.ss_starts[ss_idx_on_query]) < ALLOWED_ERROR_BASES
                && std::abs(ref.ss_ends[ss_idx_on_ref] - query.ss_ends[ss_idx_on_query]) < ALLOWED_ERROR_BASES) {
                ss_matches.insert(ss_idx_on_ref);
            }
        }
    }
    return ss_matches.size();
}
std::size_t n_overlapping_exons(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query)
{
    if (ref.is_unaligned() || query.is_unaligned() || ref.contig_name != query.contig_name) {
        return 0;
    }
    std::vector<std::size_t> hit_ids {};

    std::set<std::size_t> exon_matches;

    for (decltype(query.exon_starts.size()) ss_idx_on_query = 0; ss_idx_on_query < query.exon_starts.size();
        ++ss_idx_on_query) {
        ref.ss_with_err.overlap(query.exon_starts[ss_idx_on_query], query.exon_ends[ss_idx_on_query], hit_ids);
        for (const auto& ss_idx_on_ref : hit_ids) {
            if (std::abs(ref.exon_starts[ss_idx_on_ref] - query.exon_starts[ss_idx_on_query]) < ALLOWED_ERROR_BASES
                && std::abs(ref.exon_ends[ss_idx_on_ref] - query.exon_ends[ss_idx_on_query]) < ALLOWED_ERROR_BASES) {
                exon_matches.insert(ss_idx_on_ref);
            }
        }
    }
    return exon_matches.size();
}
} // namespace labw::art_modern
