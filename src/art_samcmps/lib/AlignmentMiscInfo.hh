#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <ewah/ewah.h>

#include <IITree.h>

#include <cstdlib>
#include <string>
#include <vector>

namespace labw::art_modern {
/**
 * Miscellaneous stuff that helps to compare two transcript alignments.
 */
class AlignmentMiscInfo {
public:
    DELETE_COPY(AlignmentMiscInfo)

    DEFAULT_MOVE(AlignmentMiscInfo)

    // Coordinates respect to the 5' end of the chromosome.
    std::string name;
    /**
     * This field must be set for aligned alignments.
     * Will be empty if not aligned.
     */
    std::string contig_name;
    hts_pos_t start = 0;
    hts_pos_t end = 0;
    std::vector<hts_pos_t> ss_starts;
    std::vector<hts_pos_t> ss_ends;
    std::vector<hts_pos_t> deletion_starts;
    std::vector<hts_pos_t> deletion_ends;
    std::vector<hts_pos_t> insertion_pos;
    std::vector<hts_pos_t> insertion_lengths;
    std::vector<hts_pos_t> exon_starts;
    std::vector<hts_pos_t> exon_ends;
    std::size_t exonic_length = 0;
    ewah::EWAHBoolArray<std::size_t> base_set; // This have to be unsigned.
    IITree<hts_pos_t, std::size_t> ss;
    IITree<hts_pos_t, std::size_t> deletions;
    IITree<hts_pos_t, std::size_t> exons;
    IITree<hts_pos_t, std::size_t> ss_with_err;
    IITree<hts_pos_t, std::size_t> deletions_with_err;
    IITree<hts_pos_t, std::size_t> exons_with_err;

    [[nodiscard]] bool is_unaligned() const;

    AlignmentMiscInfo() = default;

    ~AlignmentMiscInfo() = default;

    void from_bam_record(const bam_hdr_t* hdr, const bam1_t* aln);

    void from_gffutils_bed_line(const std::string& line);

    void reset();

private:
    void update_cr_index();
};

/**
 * Count number of overlapping bases.
 * @param ref The ground-truth alignment.
 * @param query The alignment aligned by the aligner.
 * @return As described.
 */
std::size_t n_overlapping_base(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query);

/**
 * Count number of overlapping splice sites.
 *
 * Overlaps are determined using ALLOWED_ERROR_BASES.
 *
 * @param ref The ground-truth alignment.
 * @param query The alignment aligned by the aligner.
 * @return As described.
 */
std::size_t n_overlapping_ss(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query);

/**
 * Count number of overlapping exons.
 *
 * Overlaps are determined using ALLOWED_ERROR_BASES.
 *
 * @param ref The ground-truth alignment.
 * @param query The alignment aligned by the aligner.
 * @return As described.
 */
std::size_t n_overlapping_exons(const AlignmentMiscInfo& ref, const AlignmentMiscInfo& query);
} // namespace labw::art_modern
