#include "art/ArtRead.hh"

#include "art/random_generator.hh"

#include "art_modern_config.h" // NOLINT: For CEU_IS_DEBUG
#include "libam/Constants.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/utils/seq_utils.hh"

#include <htslib/hts.h>

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <utility>

namespace labw::art_modern {

void ArtRead::generate_pairwise_aln()
{
    std::size_t pos_on_aln_str = 0;
    std::size_t pos_on_read = 0;
    std::size_t pos_on_ref = 0;
    const std::size_t maxk = art_params_.read_len + 1 + 1 + indel_.size();
    aln_read_.resize(maxk);
    aln_ref_.resize(maxk);
    std::size_t num_match = 0;

    for (const auto& [this_pos_on_aln_str, indel] : indel_) {
        num_match = this_pos_on_aln_str - pos_on_aln_str;
        std::memcpy(aln_read_.data() + pos_on_aln_str, seq_read_.data() + pos_on_read, num_match);
        std::memcpy(aln_ref_.data() + pos_on_aln_str, seq_ref_.data() + pos_on_ref, num_match);
        pos_on_read += num_match;
        pos_on_ref += num_match;
        pos_on_aln_str = this_pos_on_aln_str;
        if (indel == ALN_GAP) {
            aln_read_[pos_on_aln_str] = ALN_GAP;
            aln_ref_[pos_on_aln_str] = seq_ref_[pos_on_ref];
            pos_on_ref++;
        } else {
            aln_read_[pos_on_aln_str] = indel;
            aln_ref_[pos_on_aln_str] = ALN_GAP;
            pos_on_read++;
        }
        pos_on_aln_str++;
    }
    num_match = seq_ref_.size() - pos_on_ref;
    std::memcpy(aln_read_.data() + pos_on_aln_str, seq_read_.data() + pos_on_read, num_match);
    std::memcpy(aln_ref_.data() + pos_on_aln_str, seq_ref_.data() + pos_on_ref, num_match);
    pos_on_aln_str += num_match;
#if (0) // Old version for historical purposes
    while (pos_on_ref < seq_ref_.size()) {
        const auto find = indel_.find(pos_on_aln_str);
        if (find == indel_.end()) { // No indel
            aln_read_[pos_on_aln_str] = seq_read_[pos_on_read];
            aln_ref_[pos_on_aln_str] = seq_ref_[pos_on_ref];
            pos_on_read++;
            pos_on_ref++;
        } else if (find->second == ALN_GAP) { // Deletion
            aln_read_[pos_on_aln_str] = ALN_GAP;
            aln_ref_[pos_on_aln_str] = seq_ref_[pos_on_ref];
            pos_on_ref++;
        } else { // Insertion
            aln_read_[pos_on_aln_str] = find->second;
            aln_ref_[pos_on_aln_str] = ALN_GAP;
            pos_on_read++;
        }
        pos_on_aln_str++;
    }
    while (true) { // Insertions after reference
        const auto find = indel_.find(pos_on_aln_str);
        if (find == indel_.end()) {
            break;
        }
        aln_read_[pos_on_aln_str] = find->second;
        aln_ref_[pos_on_aln_str] = ALN_GAP;
        pos_on_aln_str++;
    }
#endif
    aln_read_.resize(pos_on_aln_str);
    aln_ref_.resize(pos_on_aln_str);
}

void ArtRead::generate_snv_on_qual(const bool is_first_read)
{
    if (!art_params_.sep_flag) {
        art_params_.qdist.get_read_qual(qual_, art_params_.read_len, rprob_, is_first_read);
    } else if (is_first_read) {
        art_params_.qdist.get_read_qual_sep_1(qual_, seq_read_, rprob_);
    } else {
        art_params_.qdist.get_read_qual_sep_2(qual_, seq_read_, rprob_);
    }
    char achar = 0;
    rprob_.r_probs();
    for (decltype(qual_.size()) i = 0; i < qual_.size(); i++) {
        if (seq_read_[i] == 'N') {
            qual_[i] = MIN_QUAL;
            continue;
        }
        if (rprob_.tmp_probs_[i] < art_params_.err_prob[qual_[i]]) {
            achar = seq_read_[i];
            while (seq_read_[i] == achar) {
                achar = rprob_.rand_base();
            }
            seq_read_[i] = achar;
        }
    }
}

int ArtRead::generate_indels(const bool is_read_1)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    int i = 0;
    int j = 0;
    int pos = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
    // deletion
    rprob_.r_probs(per_base_del_rate.size());
    for (i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= rprob_.tmp_probs_[i]) {
            del_len = i + 1;
            j = i;
            while (j >= 0) {
                pos = rprob_.rand_pos_on_read_not_head_and_tail();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }
    rprob_.r_probs(per_base_ins_rate.size());
    for (i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if (art_params_.read_len - del_len - ins_len < i + 1) {
            continue; // ensure that there are stil enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= rprob_.tmp_probs_[i]) {
            ins_len = i + 1;
            j = i;
            while (j >= 0) {
                pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }
    return ins_len - del_len;
}

int ArtRead::generate_indels_2(const bool is_read_1)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    int pos = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

    rprob_.r_probs(per_base_ins_rate.size());
    for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= rprob_.tmp_probs_[i]) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }

    // deletion
    rprob_.r_probs(per_base_del_rate.size());
    for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if (art_params_.read_len - del_len - ins_len < i + 1) {
            continue; // ensure that enough unchanged position for mutation
        }

        if (per_base_del_rate[i] >= rprob_.tmp_probs_[i]) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                pos = rprob_.rand_pos_on_read_not_head_and_tail();
                if (pos == 0) {
                    continue;
                }
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }
    return ins_len - del_len;
}

void ArtRead::ref2read(std::string seq_ref, const bool is_plus_strand, const hts_pos_t pos_on_contig)
{
    pos_on_contig_ = pos_on_contig;
    is_plus_strand_ = is_plus_strand;
    seq_ref_ = std::move(seq_ref);
    normalize_inplace(seq_ref_);
    if (!is_plus_strand) {
        revcomp_inplace(seq_ref_);
    }
    std::size_t pos_on_aln_str = 0;
    std::size_t pos_on_read = 0;
    std::size_t pos_on_ref = 0;
    const std::size_t maxk = art_params_.read_len + 1 + 1 + indel_.size();
    aln_read_.resize(maxk);
    aln_ref_.resize(maxk);
    std::size_t num_match = 0;

    for (const auto& [this_pos_on_aln_str, indel] : indel_) {
        num_match = this_pos_on_aln_str - pos_on_aln_str;
        std::memcpy(seq_read_.data() + pos_on_read, seq_ref_.data() + pos_on_ref, num_match);
        pos_on_read += num_match;
        pos_on_ref += num_match;
        pos_on_aln_str = this_pos_on_aln_str;
        if (indel == ALN_GAP) {
            pos_on_ref++;
        } else {
            seq_read_[pos_on_read] = indel;
            pos_on_read++;
        }
        pos_on_aln_str++;
    }
    num_match = seq_ref_.size() - pos_on_ref;
    std::memcpy(seq_read_.data() + pos_on_read, seq_ref_.data() + pos_on_ref, num_match);
#if (0) // Old code for historical purposes
    for (decltype(seq_ref_.size()) pos_on_ref = 0; pos_on_ref < seq_ref_.size();) {
        const auto find = indel_.find(k);
        if (find == indel_.end()) { // No indel
            seq_read_[pos_on_read] = seq_ref_[pos_on_ref];
            pos_on_ref++;
            pos_on_read++;
        } else if (find->second == ALN_GAP) { // Deletion
            pos_on_ref++;
        } else { // Insertion
            seq_read_[pos_on_read] = find->second;
            pos_on_read++;
        }
        k++;
    }
    while (true) { // Insertions after reference
        const auto find = indel_.find(k);
        if (find == indel_.end()) {
            break;
        }
        seq_read_[pos_on_read] = find->second;
        pos_on_read++;
        k++;
    }
#endif
#ifdef CEU_IS_DEBUG
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        throw ReadGenerationException("Generated seq_read_ with unequal sizes");
    }
#endif
}

ArtRead::ArtRead(
    const ArtParams& art_params, const std::string& contig_name, const std::string& read_name, Rprob& rprob)
    : art_params_(art_params)
    , contig_name_(contig_name)
    , read_name_(read_name)
    , rprob_(rprob)
{
    seq_read_.resize(art_params_.read_len);
    qual_.resize(art_params_.read_len);
}

PairwiseAlignment ArtRead::to_pwa()
{
    std::string qual_str = qual_to_str(qual_);
    return { std::move(read_name_), std::move(contig_name_), std::move(seq_read_), std::move(seq_ref_),
        std::move(qual_str), std::move(aln_read_), std::move(aln_ref_), pos_on_contig_, is_plus_strand_ };
}
bool ArtRead::is_good() const
{
    if (std::count(seq_read_.begin(), seq_read_.end(), 'N') > art_params_.max_n) {
        return false;
    }
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        goto error;
    }
    if (static_cast<int>(qual_.size()) != art_params_.read_len) {
        goto error;
    }
    return true;
error:
    return false; // FIXME: No idea why this occurs.
}

} // namespace labw::art_modern
