#pragma once

#include <htslib/hts.h>

#include <string>
#include <vector>

namespace labw::art_modern {
class Transcript {
public:
    /**
     * Splice site position on 5' to 3' of the genome.
     */
    const std::vector<hts_pos_t> splice_site_positions;
    /**
     * Splice site length on 5' to 3' of the genome.
     */
    const std::vector<hts_pos_t> splice_site_lengths;
    /**
     * Whether the transcript is on the reverse strand.
     */
    const bool is_reverse;
    /** Transcript start position in genome.
     * May be the 1st base of the 1st exon (on positive strand)
     * or the last base on the last exon (on negative strand).
     * Is 0-based inclusive.
     */
    const hts_pos_t start;
    const hts_pos_t end;
    /**
     * TID of the contig on the reference genome FASTA.
     */
    const int tid_on_genome;
    /**
     * TID of the transcript on the reference transcriptome FASTA.
     */
    const int tid_on_transcriptome;
    /**
     * Transcript ID.
     */
    const std::string transcript_id;
    const std::string contig_name;
    const hts_pos_t unspliced_length;

    static Transcript from_gffread_bed_line(const std::string& line, sam_hdr_t* thdr, sam_hdr_t* ghdr);
};
} // namespace labw::art_modern
