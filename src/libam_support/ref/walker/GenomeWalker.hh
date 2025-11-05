#pragma once
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <htslib/hts.h>

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>

namespace labw::art_modern {

enum class VariantType : std::uint8_t { SNV, INS, DEL, BND };

/**
 * A walker to walk along the genome.
 */
class GenomeWalker {
public:
    explicit GenomeWalker(std::shared_ptr<BaseFastaFetch>&& fasta_fetcher);
    explicit GenomeWalker(const std::shared_ptr<BaseFastaFetch>& fasta_fetcher);

    void attach(std::size_t seq_id, hts_pos_t start, bool on_pos_strand = true);

    /**
     * Where the walker is currently at.
     */
    [[nodiscard]] std::size_t get_seq_id() const;
    /**
     * Where the walker is currently at.
     */
    [[nodiscard]] hts_pos_t get_position() const;
    /**
     * Where the walker is currently at.
     */
    [[nodiscard]] bool is_on_pos_strand() const;

    /**
     * Reset the walker.
     */
    void clear();

    /**
     * Skip the given offset, typically due to an intron. Variants inside the intron will be considered.
     *
     * @param offset As described.
     */
    void skip(hts_pos_t offset);

    /**
     * Move forward the given offset without considering variants.
     *
     * @param offset As described.
     */
    void forward(hts_pos_t offset);

    /**
     * Walk the given length of bases. Variants will be considered.
     *
     * @param length As described.
     * @return The bases.
     */
    std::string walk(hts_pos_t length);

    /**
     * Check whether the walker has reached the genome boundary.
     *
     * @return As described.
     */
    [[nodiscard]] bool reached_genome_boundary() const;

private:
    constexpr static hts_pos_t VARIANT_CACHE_SIZE_ = 4096; // Number of bases to cache variants for.
    void* variant_cache_; // Placeholder for variant cache, implementation not shown here.

    /**
     * Refresh the variant cache of next VARIANT_CACHE_SIZE_ bases.
     */
    void refresh_variant_cache_();

    std::shared_ptr<BaseFastaFetch> fasta_fetcher_;
    std::size_t seq_id_ = 0;
    hts_pos_t position_ = 0;
    bool on_pos_strand_ = true;
};

} // namespace labw::art_modern
