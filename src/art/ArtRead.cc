#include "ArtRead.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"

using namespace std;

namespace labw::art_modern {

void ArtRead::generate_pairwise_aln()
{
    int k = 0;
    int pos_on_read = 0;
    std::ostringstream aln_seq_ss;
    std::ostringstream aln_ref_ss;
    for (decltype(seq_ref.size()) pos_on_ref = 0; pos_on_ref < seq_ref.size();) {
        if (indel_.find(k) == indel_.end()) { // No indel
            aln_seq_ss << seq_read_[pos_on_read];
            aln_ref_ss << seq_ref[pos_on_ref];
            pos_on_read++;
            pos_on_ref++;
        } else if (indel_[k] == ALN_GAP) { // Deletion
            aln_seq_ss << ALN_GAP;
            aln_ref_ss << seq_ref[pos_on_ref];
            pos_on_ref++;
        } else { // Insertion
            aln_seq_ss << indel_[k];
            aln_ref_ss << ALN_GAP;
            pos_on_read++;
        }
        k++;
    }
    while (indel_.find(k) != indel_.end()) { // Insertions after reference
        aln_seq_ss << indel_[k];
        aln_ref_ss << ALN_GAP;
        pos_on_read++;
        k++;
    }
    aln_read_ = aln_seq_ss.str();
    aln_ref_ = aln_ref_ss.str();
}

void ArtRead::generate_snv_on_qual(bool is_first_read)
{
    if (!art_params_.sep_flag) {
        art_params_.qdist.get_read_qual(qual_, art_params_.read_len, rprob_, is_first_read);
    } else if (is_first_read) {
        art_params_.qdist.get_read_qual_sep_1(qual_, seq_read_, rprob_);
    } else {
        art_params_.qdist.get_read_qual_sep_2(qual_, seq_read_, rprob_);
    }
    char achar;
    const auto probs = rprob_.r_probs(qual_.size());
    for (decltype(qual_.size()) i = 0; i < qual_.size(); i++) {
        if (seq_read_[i] == 'N') {
            qual_[i] = MIN_QUAL;
            continue;
        }
        if (probs[i] < art_params_.err_prob[qual_[i]]) {
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
    int i;
    int j;
    int pos;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
    // deletion
    const auto probs_del = rprob_.r_probs(per_base_ins_rate.size());
    for (i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= probs_del[i]) {
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
    const auto probs_ins = rprob_.r_probs(per_base_ins_rate.size());
    for (i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= probs_ins[i]) {
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
    return (ins_len - del_len);
};

int ArtRead::generate_indels_2(const bool is_read_1)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

    const auto probs_ins = rprob_.r_probs(per_base_ins_rate.size());
    for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= probs_ins[i]) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }

    // deletion
    const auto probs_del = rprob_.r_probs(per_base_ins_rate.size());
    for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }

        if (per_base_del_rate[i] >= probs_del[i]) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
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
    return (ins_len - del_len);
};

void ArtRead::ref2read()
{
    int k = 0;
    int pos_on_read = 0;
    for (decltype(seq_ref.size()) pos_on_ref = 0; pos_on_ref < seq_ref.size();) {
        if (indel_.find(k) == indel_.end()) { // No indel
            seq_read_[pos_on_read] = seq_ref[pos_on_ref];
            pos_on_ref++;
            pos_on_read++;
        } else if (indel_[k] == ALN_GAP) { // Deletion
            pos_on_ref++;
        } else { // Insertion
            seq_read_[pos_on_read] = indel_[k];
            pos_on_read++;
        }
        k++;
    }
    while (indel_.find(k) != indel_.end()) { // Insertions after reference
        seq_read_[pos_on_read] = indel_[k];
        pos_on_read++;
        k++;
    }
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        throw ReadGenerationException("Generated seq_read_ with unequal sizes");
    }
}

ArtRead::ArtRead(
    const ArtParams& art_params, Rprob& rprob, const std::string& contig_name, const std::string& read_name)
    : art_params_(art_params)
    , rprob_(rprob)
    , read_name_(read_name)
    , contig_name_(contig_name)
{
    seq_read_.resize(art_params_.read_len);
    seq_ref.reserve(art_params_.read_len);
    qual_.resize(art_params_.read_len);
}

PairwiseAlignment ArtRead::to_pwa() const
{
    return { read_name_, contig_name_, seq_read_, seq_ref, qual_to_str(qual_), aln_read_, aln_ref_, pos_on_contig,
        is_plus_strand };
}
bool ArtRead::is_good() const
{
    if (std::count(seq_read_.begin(), seq_read_.end(), 'N') > 0) { // TODO: Add params back.
        return false;
    }
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        goto error;
    }
    if (static_cast<int>(qual_.size()) != art_params_.read_len) {
        goto error;
    }
    if (aln_read_.size() != aln_ref_.size()) {
        goto error;
    }
    return true;
error:
    // #ifdef CEU_CM_IS_DEBUG
    //         abort_mpi();
    // #else
    return false; // FIXME: No idea why this occurs.
    // #endif
}

} // namespace labw::art_modern; // namespace labw