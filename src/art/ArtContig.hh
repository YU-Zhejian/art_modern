#pragma once

#include <set>
#include <string>

#include "ArtParams.hh"
#include "ArtRead.hh"
#include "fasta/BaseFastaFetch.hh"
#include "random_generator.hh"

namespace labw::art_modern {

class ArtContig {

public:
    ArtContig(BaseFastaFetch* fasta_fetch, size_t seq_id, const ArtParams& art_params, Rprob& rprob);

    void generate_read_se(bool is_plus_strand, ArtRead& read_1);

    void generate_read_pe(bool is_plus_strand, bool is_mp, ArtRead& read_1, ArtRead& read_2);

    const std::string seq_name;

private:
    hts_pos_t generate_fragment_length() const;
    BaseFastaFetch* fasta_fetch_;
    const ArtParams& art_params_;
    Rprob& rprob_;
    const size_t seq_id_;
    const hts_pos_t valid_region_;
    const hts_pos_t ref_len_;
};

} // namespace labw::art_modern // namespace labw