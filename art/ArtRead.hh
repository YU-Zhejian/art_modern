#pragma once

#include <map>
#include <string>

#include "ArtParams.hh"
#include "Empdist.hh"
#include "PairwiseAlignment.hh"
#include "random_generator.hh"
#include "seq_utils.hh"

namespace labw {
namespace art_modern {

    struct TooMuchNException : public std::exception {
        const char* what() const noexcept override { return "Too much N in the contig"; }
    };

    class ArtRead {
    public:
        ArtRead(const ArtParams& art_params, Rprob& rprob);

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
        void assess_num_n();

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
        const ArtParams& art_params_;
        Rprob& rprob_;
    };

    struct ArtReadPair {
        ArtRead read_1;
        ArtRead read_2;
    };

} // namespace art_modern
} // namespace labw