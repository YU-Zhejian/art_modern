#include "ArtRead.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"

using namespace std;

namespace labw {
namespace art_modern {

    void ArtRead::generate_pairwise_aln()
    {
        aln_read = seq_read;
        aln_ref = seq_ref;
        for (auto pos = 0; pos < art_params_.read_len; pos++) {
            if (indel_[pos] == 0) {
                continue;
            }
            (indel_[pos] != ALN_GAP ? aln_ref : aln_read).insert(indel_[pos], 1, ALN_GAP);
        }
    }

    void ArtRead::generate_snv_on_qual(std::vector<int>& qual)
    {
        for (auto i = 0; i < qual.size(); i++) {
            if (seq_read[i] == 'N') {
                qual[i] = 1;
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
        int ins_len = 0;
        int del_len = 0;
        std::fill(indel_.begin(), indel_.end(), 0);
        auto per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
        auto per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
        // deletion
        for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
            if (per_base_del_rate[i] >= rprob_.r_prob()) {
                del_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
                    if (indel_[pos] == 0) {
                        indel_[pos] = ALN_GAP;
                        j--;
                    }
                }
                break;
            }
        }

        for (auto i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
            if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
                continue; // ensure that enough unchanged position for mutation
            }
            if (per_base_ins_rate[i] >= rprob_.r_prob()) {
                ins_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read();
                    if (indel_[pos] == 0) {
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
        int ins_len = 0;
        int del_len = 0;
        indel_.clear();
        auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
        auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

        for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
            if (per_base_ins_rate[i] >= rprob_.r_prob()) {
                ins_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read();
                    if (indel_[pos] == 0) {
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
                    if (indel_[pos] == 0) {
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
        for (auto i = 0; i < seq_ref.size();) {
            if (indel_[k] == 0) {
                seq_read.push_back(seq_ref[i]);
                i++;
                k++;
            } else if (indel_[k] == ALN_GAP) {
                i++;
                k++;
            } else {
                seq_read.push_back(indel_[k]);
                k++;
            }
        }
        while (indel_[k] != 0) {
            seq_read.push_back(indel_[k]);
            k++;
        }
    }

    ArtRead::ArtRead(const ArtParams& art_params, Rprob& rprob)
        : art_params_(art_params)
        , rprob_(rprob)
        , indel_(art_params.read_len)
    {
    }

    void ArtRead::assess_num_n() const
    {
        auto num_n = std::count(seq_read.begin(), seq_read.end(), 'N');

        if (num_n > 0) { // TODO: Add params back.
            throw TooMuchNException();
        }
    }

} // namespace art_modern
}; // namespace labw