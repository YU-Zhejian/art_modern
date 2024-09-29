#pragma once

#include <set>
#include <string>

#include "ArtParams.hh"
#include "ArtRead.hh"
#include "fasta/BaseFastaFetch.hh"
#include "random_generator.hh"

namespace labw {
namespace art_modern {

    class ArtContig {

    public:
        ArtContig(BaseFastaFetch* fasta_fetch, size_t seq_id, const ArtParams& art_params, Rprob& rprob);

        ArtRead generate_read_se(bool is_plus_strand);

        ArtReadPair generate_read_pe(bool is_plus_strand);

        ArtReadPair generate_read_mp(bool is_plus_strand);

        BaseFastaFetch* fasta_fetch_;
        const ArtParams& art_params_;
        Rprob& rprob_;
        const size_t seq_id_;
        const hts_pos_t valid_region_;
        const hts_pos_t ref_len_;
        const std::string id_;
    };

} // namespace art_modern
} // namespace labw