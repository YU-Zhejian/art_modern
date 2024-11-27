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
    bool is_plus_strand = false;
    long pos_on_contig;
    std::string seq_ref;

    // Disable copy constructors
    ArtRead(const ArtRead&) = delete;
    ArtRead& operator=(ArtRead&&) = delete;
    // Enable move constructors
    ArtRead(ArtRead&& other) noexcept = default;

    ArtRead(const ArtParams& art_params, Rprob& rprob, const std::string& contig_name, const std::string& read_name);
    PairwiseAlignment to_pwa() const;

    int generate_indels(bool is_read_1);
    // number of deletions <= number of insertions
    int generate_indels_2(bool is_read_1);

    /**
     * Populate the read while adding insertions and deletions.
     */
    void ref2read();

    void generate_pairwise_aln();

    /**
     * Add point mutations to random bases based on empirical dist of quali scores
     */
    void generate_snv_on_qual(bool is_first_read);
    bool is_good() const;

private:
    const ArtParams& art_params_;
    Rprob& rprob_;
    std::map<int, char, std::less<>> indel_;
    const std::string read_name_;
    const std::string contig_name_;
    std::string aln_read_;
    std::string aln_ref_;
    std::vector<int> qual_;
    std::string seq_read_;
};

} // namespace labw::art_modern // namespace labw