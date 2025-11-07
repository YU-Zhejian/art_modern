#include "tsam2gsam/lib/Transcript.hh"

#include "tsam2gsam/lib/TidNotFound.hh"
#include "tsam2gsam/lib/cyh_proj_utils.hh"

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include <htslib/sam.h>

#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
Transcript Transcript::from_gffread_bed_line(const std::string& line, sam_hdr_t* thdr, sam_hdr_t* ghdr)
{
    std::string transcript_id;
    std::string contig_name;
    std::vector<int32_t> splice_site_positions;
    splice_site_positions.reserve(NUM_EXONS_EXPECTED);
    std::vector<int32_t> splice_site_lengths;
    splice_site_lengths.reserve(NUM_EXONS_EXPECTED);
    bool is_reverse = false;
    int32_t t_start = 0;
    int32_t t_end = 0;
    int tid_on_transcriptome = 0;
    int tid_on_genome = 0;
    std::array<std::string, NUM_GFFUTILS_BED_FIELDS> tokens;
    std::string token;
    int n_exons = 0;
    int32_t accumulated_exon_length = 0;
    int32_t accumulated_intron_length = 0;
    int32_t current_exon_length = 0;
    int32_t current_intron_length = 0;
    int32_t current_exon_pos_on_unspliced_transcript = 0;

    std::istringstream line_ss(line);
    for (int i = 0; i < NUM_GFFUTILS_BED_FIELDS; ++i) {
        std::getline(line_ss, token, '\t');
        tokens[i] = std::move(token);
    }
    contig_name = std::move(tokens[0]);
    t_start = std::stoi(tokens[1]);
    t_end = std::stoi(tokens[2]);
    transcript_id = std::move(tokens[3]);
    is_reverse = tokens[5] == "-";

    n_exons = std::stoi(tokens[9]);
    std::istringstream elen_ss(tokens[10]);
    std::istringstream epos_ss(tokens[11]);
    accumulated_exon_length = 0;
    accumulated_intron_length = 0;
    for (auto exon_id = 0; exon_id < n_exons; ++exon_id) {
        std::getline(epos_ss, token, ',');
        current_exon_pos_on_unspliced_transcript = std::stoi(token);
        std::getline(elen_ss, token, ',');
        current_exon_length = std::stoi(token);

        current_intron_length
            = current_exon_pos_on_unspliced_transcript - accumulated_exon_length - accumulated_intron_length;
        if (current_intron_length != 0) {
            splice_site_lengths.emplace_back(current_intron_length);
            splice_site_positions.emplace_back(accumulated_exon_length);
        }
        accumulated_exon_length += current_exon_length;
        accumulated_intron_length += current_intron_length;
    }
    if (thdr == nullptr) {
        tid_on_transcriptome = -1;
    } else {
        tid_on_transcriptome = sam_hdr_name2tid(thdr, transcript_id.c_str());
        if (tid_on_transcriptome < 0) {
            throw TidNotFound();
        }
    }
    if (ghdr == nullptr) {
        tid_on_genome = -1;
    } else {
        tid_on_genome = sam_hdr_name2tid(ghdr, contig_name.c_str());
        if (tid_on_genome < 0) {
            throw TidNotFound();
        }
    }
#ifdef CEU_CM_IS_DEBUG
    if (t_end - t_start != accumulated_exon_length + accumulated_intron_length) {
        std::abort();
    }
#endif
    return { splice_site_positions, splice_site_lengths, is_reverse, t_start, t_end, tid_on_genome,
        tid_on_transcriptome, transcript_id, contig_name, accumulated_exon_length };
}
} // namespace labw::art_modern
