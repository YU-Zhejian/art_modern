#pragma once

#include <set>
#include <string>

#include "ArtParams.hh"
#include "ArtRead.hh"

namespace labw {
namespace art_modern {

    struct ArtReadPair {
        ArtRead read_1;
        ArtRead read_2;
    };

    class ArtContig {

    public:
        ArtContig(std::string ref_seq, std::string id, const ArtParams& art_params);

        ArtRead generate_read_se() const;

        ArtReadPair generate_read_pe() const;

        ArtReadPair generate_read_mp() const;

        void mask_n_region(int max_num_n);

        ArtParams _art_params;
        /**
         * Current reference genome contig.
         * Contains upper-case characters only.
         */
        std::string _ref_seq;
        std::string _id;
        /**
         * Reverse-complementary of {@link _ref_seq}.
         */
        std::string _ref_seq_cmp;
        int _valid_region;
        /**
         * Region with n
         */
        std::set<size_t> _masked_pos;
    };

} // namespace art_modern
} // namespace labw