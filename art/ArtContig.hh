#pragma once

#include <set>
#include <string>

#include "ArtParams.hh"
#include "ArtRead.hh"

namespace labw {
namespace art_modern {

    class ArtContig {

    public:
        ArtContig(std::string ref_seq, std::string id, const ArtParams& art_params);

        ArtRead generate_read_se() const;

        ArtReadPair generate_read_pe() const;

        ArtReadPair generate_read_mp() const;

        ArtParams _art_params;
        /**
         * Current reference genome contig.
         * Contains upper-case characters only.
         */
        std::string _ref_seq;
        std::string _id;
        int _valid_region;
    };

} // namespace art_modern
} // namespace labw