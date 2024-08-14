#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ArtParams.hh"
#include "Empdist.hh"
#include "PairwiseAlignment.hh"
#include "seq_utils.hh"

namespace labw {
namespace art_modern {

    class ArtRead {
    public:
        explicit ArtRead(ArtParams art_params);

        std::map<int, char> indel;
        std::map<int, char> substitution;
        bool is_plus_strand = false;
        long bpos;
        std::string seq_read;
        std::string seq_ref;
        std::string aln_read;
        std::string aln_ref;

        int generate_indels(int read_len, bool is_read_1);
        // number of deletions <= number of insertions
        int generate_indels_2(int read_len, bool is_read_1);

        /**
         * Populate the read while adding insertions and deletions.
         */
        void ref2read();

        void generate_pairwise_aln();

        /**
         * Add point mutations to random bases based on empirical dist of quali scores
         * @param qual As described.
         */
        std::vector<int> generate_snv_on_qual(const std::vector<int>& qual);

    private:
        ArtParams _art_params;
    };

    struct ArtReadPair {
        ArtRead read_1;
        ArtRead read_2;
    };

} // namespace art_modern
} // namespace labw