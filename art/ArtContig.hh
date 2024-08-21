#pragma once

#include <set>
#include <string>

#include "ArtParams.hh"
#include "ArtRead.hh"
#include "fasta/BaseFastaFetch.hh"
#include "random_generator.hh"
#include <htslib/sam.h>
#include <memory>

namespace labw {
namespace art_modern {

    class ArtContig {

    public:
        ArtContig(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, std::string id, const ArtParams& art_params,
            Rprob& rprob);

        ArtRead generate_read_se(bool is_plus_strand);

        ArtReadPair generate_read_pe(bool is_plus_strand);

        ArtReadPair generate_read_mp(bool is_plus_strand);

        const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
        const ArtParams& art_params_;
        Rprob& rprob_;
        const std::string id_;
        const hts_pos_t valid_region_;
        const hts_pos_t ref_len_;
    };

} // namespace art_modern
} // namespace labw