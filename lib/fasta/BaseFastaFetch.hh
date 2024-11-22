#pragma once
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <string>
#include <vector>

namespace labw::art_modern {
/**
 * Random access FASTA file.
 */
class BaseFastaFetch {
public:
    BaseFastaFetch(const BaseFastaFetch&) = delete;
    BaseFastaFetch(BaseFastaFetch&&) = delete;
    BaseFastaFetch& operator=(const BaseFastaFetch&) = delete;
    BaseFastaFetch& operator=(BaseFastaFetch&&) = delete;
    BaseFastaFetch();

    explicit BaseFastaFetch(const std::tuple<std::vector<std::string>, std::vector<hts_pos_t>>& seq_names_lengths);
    BaseFastaFetch(std::vector<std::string> seq_names, std::vector<hts_pos_t> seq_lengths);

    /**
     * This method is thread-safe since mutex is used for non-thread-safe implementations.
     *
     * @param seq_id Contig name.
     * @param start 0-based inclusive start point.
     * @param end 0-based exclusive end point.
     * @return Fetched sequence.
     */
    virtual std::string fetch(std::size_t seq_id, hts_pos_t start, hts_pos_t end) = 0;

    /**
     * This method is thread-safe since mutex is used for non-thread-safe implementations.
     * Fetch the entire contig.
     *
     * @param seq_id Contig name.
     * @return Fetched sequence.
     */
    virtual std::string fetch(std::size_t seq_id);

    /**
     * Default destructor.
     */
    virtual ~BaseFastaFetch();

    void update_sam_header(sam_hdr_t* header) const;

    /**
     * Get the length of the desired sequence.
     * This operation should take $O(1)$ complexity.
     *
     * @param seq_id As described.
     * @return As described.
     */
    hts_pos_t seq_len(std::size_t seq_id) const;

    /**
     * Get the name of the desired sequence.
     * @param seq_id As described.
     * @return As described.
     */
    std::string seq_name(std::size_t seq_id) const;

    /**
     * Get number of sequences inside.
     *
     * @return As described.
     */
    size_t num_seqs() const;

    bool empty() const;

protected:
    const std::vector<std::string> seq_names_;
    const std::vector<hts_pos_t> seq_lengths_;
};
}
