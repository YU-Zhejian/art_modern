#pragma once

#include "art/ArtParams.hh"
#include "art/random_generator.hh"

#include "libam/Dtypes.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/utils/class_macros_utils.hh"

#include <absl/base/attributes.h>

#include <htslib/hts.h>

#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace labw::art_modern {

struct ReadGenerationException : public std::runtime_error {
    using runtime_error::runtime_error;
} ABSL_ATTRIBUTE_PACKED;

class ArtRead {
public:
    // Disable constructors
    DELETE_MOVE(ArtRead)
    DELETE_COPY(ArtRead)
    ~ArtRead() = default;

    ArtRead(const ArtParams& art_params, std::string contig_name, const std::string& read_name, Rprob& rprob);
    [[nodiscard]] PairwiseAlignment to_pwa();

    int generate_indels(bool is_read_1);
    // number of deletions <= number of insertions
    int generate_indels_2(bool is_read_1);

    /*!
     * Populate the read while adding insertions and deletions.
     */
    void ref2read(std::string seq_ref, bool is_plus_strand, hts_pos_t pos_on_contig);

    /*!
     * Populate aln_read_ and aln_ref_
     *
     * TODO: This stuff is not needed for FASTQ. Need to modify the control flow.
     */
    void generate_pairwise_aln();

    /*!
     * Add point mutations to random bases based on empirical dist of quali scores
     */
    void generate_snv_on_qual(bool is_first_read);
    [[nodiscard]] bool is_good() const;

private:
    std::string aln_read_;
    std::string aln_ref_;
    const ArtParams& art_params_;
    std::string contig_name_;
    std::map<int, char, std::less<>> indel_;
    bool is_plus_strand_ = false;
    hts_pos_t pos_on_contig_ = 0;
    std::vector<am_qual_t> qual_;
    std::string read_name_;
    Rprob& rprob_;
    std::string seq_read_;
    std::string seq_ref_;
};

} // namespace labw::art_modern