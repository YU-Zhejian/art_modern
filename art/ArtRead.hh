#pragma once

#include <map>
#include <string>
#include <vector>

#include "ArtParams.hh"
#include "Empdist.hh"
#include "PairwiseAlignment.hh"
#include "seq_utils.hh"

namespace labw {
namespace art_modern {

    class ArtCustomException : std::exception {
    public:
        explicit ArtCustomException(const std::string& cause)
            : _cause(cause) {};

        const char* what() const noexcept override
        {
            return _cause.c_str();
        }

    private:
        std::string _cause;
    };

    class ArtRead {
    public:
        explicit ArtRead(const ArtParams& art_params)
            : _art_params(art_params)
        {
        }

        std::map<int, char> indel;
        std::map<int, char> substitution;
        bool is_plus_strand = false;
        long bpos = 0; // parent
        std::string seq_read;
        std::string seq_ref;

        int generate_indels(int read_len, bool is_read_1);
        // number of deletions <= number of insertions
        int generate_indels_2(int read_len, bool is_read_1);

        /**
         * Populate the read while adding insertions and deletions.
         */
        void ref2read();

        PairwiseAlignment generate_pairwise_aln() const;

        /**
         * Add point mutations to random bases based on empirical dist of quali scores
         * @param qual As described.
         */
        std::vector<int> generate_snv_on_qual(const std::vector<int>& qual);

    private:
        ArtParams _art_params;
    };

}
}