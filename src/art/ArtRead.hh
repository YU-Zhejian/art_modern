#pragma once

#include <map>
#include <string>

#include "ArtParams.hh"
#include "Empdist.hh"
#include "PairwiseAlignment.hh"
#include "random_generator.hh"
#include "utils/seq_utils.hh"

namespace labw::art_modern {

struct ReadGenerationException : public std::runtime_error {
    using runtime_error::runtime_error;
};

class ArtRead {
public:
    // Disable constructors
    ArtRead(const ArtRead&) = delete;
    ArtRead operator=(const ArtRead&) = delete;
    ArtRead& operator=(ArtRead&&) = delete;
    ArtRead(ArtRead&& other) noexcept = delete;

    ArtRead(const ArtParams& art_params, const std::string& contig_name, const std::string& read_name, Rprob& rprob);
    [[nodiscard]] PairwiseAlignment to_pwa();

    int generate_indels(bool is_read_1, std::vector<double>& probs_indel);
    // number of deletions <= number of insertions
    int generate_indels_2(bool is_read_1, std::vector<double>& probs_indel);

    /**
     * Populate the read while adding insertions and deletions.
     */
    void ref2read(std::string seq_ref, bool is_plus_strand, hts_pos_t pos_on_contig);

    void generate_pairwise_aln();

    /**
     * Add point mutations to random bases based on empirical dist of quali scores
     */
    void generate_snv_on_qual(bool is_first_read, std::vector<double>& tmp_qual_probs);
    [[nodiscard]] bool is_good() const;

private:
    std::string aln_read_;
    std::string aln_ref_;
    const ArtParams& art_params_;
    std::string contig_name_;
    std::map<int, char, std::less<>> indel_;
    bool is_plus_strand_ = false;
    long pos_on_contig_;
    std::vector<int> qual_;
    std::string read_name_;
    Rprob& rprob_;
    std::string seq_read_;
    std::string seq_ref_;
};

} // namespace labw::art_modern // namespace labw