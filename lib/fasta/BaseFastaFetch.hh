#pragma once
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace labw {
namespace art_modern {
    /**
     * Random access FASTA file.
     */
    class BaseFastaFetch {
    public:
        BaseFastaFetch(const BaseFastaFetch&) = delete;
        BaseFastaFetch(BaseFastaFetch&&) = delete;
        BaseFastaFetch& operator=(const BaseFastaFetch&) = delete;
        BaseFastaFetch& operator=(BaseFastaFetch&&) = delete;

        explicit BaseFastaFetch(std::unordered_map<std::string, hts_pos_t> seq_lengths);
        /**
         * This method is thread-safe since mutex is used for non-thread-safe implementations.
         *
         * @param seq_name Contig name.
         * @param start 0-based inclusive start point.
         * @param end 0-based exclusive end point.
         * @return Fetched sequence.
         */
        virtual std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end) = 0;

        /**
         * Default destructor.
         */
        virtual ~BaseFastaFetch();

        void update_sam_header(sam_hdr_t* header) const;

        /**
         * Get length of desired sequence.
         * @param seq_name As described.
         * @return As described.
         */
        hts_pos_t seq_len(const std::string& seq_name) const;

        /**
         * Get number of sequences inside.
         *
         * @return As described.
         */
        size_t num_seqs() const;
        const std::vector<std::string>& seq_names() const;

        const std::unordered_map<std::string, hts_pos_t> seq_lengths_;
        const std::vector<std::string> seq_names_;
    };
}
}
