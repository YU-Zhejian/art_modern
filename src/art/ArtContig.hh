
#pragma once

#include <cstddef>
#include <string>

#include <htslib/hts.h>

#include "art/ArtParams.hh"
#include "art/ArtRead.hh"
#include "art/random_generator.hh"

#include "libam/ref/fetch/BaseFastaFetch.hh"

namespace labw::art_modern {

class ArtContig {

public:
    ArtContig(BaseFastaFetch* fasta_fetch, size_t seq_id, const ArtParams& art_params, Rprob& rprob);

    void generate_read_se(bool is_plus_strand, ArtRead& read_1);
    void generate_read_pe(bool is_plus_strand, bool is_mp, ArtRead& read_1, ArtRead& read_2);

    /** Name of the contig */
    const std::string seq_name;
    /**Size of the contig */
    const hts_pos_t seq_size;

private:
    [[nodiscard]] hts_pos_t generate_fragment_length() const;
    const ArtParams& art_params_;
    BaseFastaFetch* fasta_fetch_;
    Rprob& rprob_;
    const size_t seq_id_;
    const hts_pos_t valid_region_;
};

} // namespace labw::art_modern