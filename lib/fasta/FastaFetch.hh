#pragma once
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <map>
#include <string>

namespace labw {
namespace art_modern {
    /**
     * Random access FASTA file.
     */
    class FastaFetch {
    public:
        /**
         * This method is thread-safe since mutex is used for non-thread-safe implementations.
         *
         * @param seq_name Contig name.
         * @param start 0-based inclusive start point.
         * @param end 0-based exclusive end point.
         * @return Fetched sequence.
         */
        virtual std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end);

        /**
         * FastaFetch::fetch with C-style API.
         *
         * @param seq_name As described.
         * @param start As described.
         * @param end As described.
         * @return The returned string should be manually freed.
         */
        virtual char* cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end);

        /**
         * Default destructor.
         */
        virtual ~FastaFetch();

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

        std::map<std::string, hts_pos_t, std::less<>> seq_lengths_;
    };
}
}
