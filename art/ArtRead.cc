#include "ArtRead.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"

using namespace std;

namespace labw::art_modern {

void ArtRead::generate_pairwise_aln()
{
    aln_read = seq_read;
    aln_ref = seq_ref;
    for (auto pos = 0; pos < art_params_.read_len; pos++) {
        if (indel_.find(pos) == indel_.end()) {
            continue;
        }
        (indel_[pos] != ALN_GAP ? aln_ref : aln_read).insert(indel_[pos], 1, ALN_GAP);
    }
}

void ArtRead::generate_snv_on_qual(std::vector<int>& qual)
{
    for (auto i = 0; i < qual.size(); i++) {
        if (seq_read[i] == 'N') {
            qual[i] = MIN_QUAL;
            continue;
        }
        if (rprob_.r_prob() < art_params_.err_prob[qual[i]]) {
            char achar = seq_read[i];
            while (seq_read[i] == achar) {
                achar = rprob_.rand_base();
            }
            seq_read[i] = achar;
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
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
    // deletion
    for (i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= rprob_.r_prob()) {
            del_len = i + 1;
            j = i;
            while (j >= 0) {
                auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }

    for (i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= rprob_.r_prob()) {
            ins_len = i + 1;
            j = i;
            while (j >= 0) {
                auto pos = rprob_.rand_pos_on_read();
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

    for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= rprob_.r_prob()) {
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
    for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }

        if (per_base_del_rate[i] >= rprob_.r_prob()) {
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
    for (auto pos_on_ref = 0; pos_on_ref < seq_ref.size();) {
        if (indel_.find(k) == indel_.end()) { // No indel
            seq_read[pos_on_read] = (seq_ref[pos_on_ref]);
            pos_on_ref++;
            pos_on_read++;
        } else if (indel_[k] == ALN_GAP) { // Deletion
            pos_on_ref++;
        } else { // Insertion
            seq_read[pos_on_read] = (indel_[k]);
            pos_on_read++;
        }
        k++;
    }
    while (indel_.find(k) != indel_.end()) { // Insertions after reference
        seq_read[pos_on_read] = (indel_[k]);
        pos_on_read++;
        k++;
    }
}

ArtRead::ArtRead(const ArtParams& art_params, Rprob& rprob)
    : art_params_(art_params)
    , rprob_(rprob)
{
    seq_read.resize(art_params_.read_len);
    seq_ref.reserve(art_params_.read_len);
}

void ArtRead::assess_num_n() const
{
    auto num_n = std::count(seq_read.begin(), seq_read.end(), 'N');

    if (num_n > 0) { // TODO: Add params back.
        throw TooMuchNException();
    }
}

} // namespace labw::art_modern; // namespace labw