#pragma once

#include <htslib/sam.h>
#include <string>
#include <vector>

namespace labw {
namespace art_modern {

    class PWAException : private std::exception {
    public:
        PWAException(std::string aligned_query, std::string aligned_ref);

        const char* what() const noexcept override;

    private:
        std::string _aligned_query;
        std::string _aligned_ref;
    };

    class PairwiseAlignment {
    public:
        PairwiseAlignment(PairwiseAlignment&& other) = delete;
        PairwiseAlignment(const PairwiseAlignment&) = delete;
        PairwiseAlignment& operator=(PairwiseAlignment&&) = delete;
        PairwiseAlignment& operator=(const PairwiseAlignment&) = delete;

        PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
            std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t align_contig_start,
            bool is_plus_strand);
        std::vector<uint32_t> generate_cigar_array(bool use_m) const;

        /**
         * Gap inserted using -
         */
        const std::string aligned_query;
        const std::string aligned_ref;
        const std::string query;
        const std::string ref;
        const std::string qual;
        const std::string read_name;
        const std::string contig_name;
        const hts_pos_t align_contig_start;
        const bool is_plus_strand;
    };

} // namespace art_modern
} // namespace labw