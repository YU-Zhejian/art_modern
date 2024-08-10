#pragma once

#include <string>

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
        PairwiseAlignment(const std::string& aligned_query,
            const std::string& aligned_ref);

        std::string generate_cigar(bool is_reverse, bool use_m) const;

    private:
        /**
         * Gap inserted using -
         */
        std::string _aligned_query;
        std::string _aligned_ref;
    };

} // namespace art_modern
} // namespace labw