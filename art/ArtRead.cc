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
        for (const auto& it : indel) {
            (it.second != ALN_GAP ? aln_ref : aln_read).insert(it.first, 1, ALN_GAP);
        }
    }

    std::vector<int> ArtRead::generate_snv_on_qual(const std::vector<int>& qual)
    {
        auto qual_mutated = qual;
        for (auto i = 0; i < qual.size(); i++) {
            if (seq_read[i] == 'N') {
                qual_mutated[i] = 1;
                continue;
            }
            if (rprob_.r_prob() < art_params_.err_prob[qual_mutated[i]]) {
                char achar = seq_read[i];
                while (seq_read[i] == achar) {
                    achar = rprob_.rand_base();
                }
                seq_read[i] = achar;
                substitution[i] = achar;
            }
        }
        return qual_mutated;
    }

    int ArtRead::generate_indels(const bool is_read_1)
    {
        int ins_len = 0;
        int del_len = 0;
        indel.clear();
        auto per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
        auto per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
        // deletion
        for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
            if (per_base_del_rate[i] >= rprob_.r_prob()) {
                del_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
                    if (indel.count(pos) == 0) {
                        indel[pos] = ALN_GAP;
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
                    if (indel.count(pos) == 0) {
                        indel[pos] = rprob_.rand_base();
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
        indel.clear();
        auto per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
        auto per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

        for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
            if (per_base_ins_rate[i] >= rprob_.r_prob()) {
                ins_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read();
                    if (indel.count(pos) == 0) {
                        indel[pos] = rprob_.rand_base();
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

            if ((art_params_.read_len - del_len - ins_len) < (i + 1))
                continue; // ensure that enough unchanged position for mutation

            if (per_base_del_rate[i] >= rprob_.r_prob()) {
                del_len = i + 1;
                for (int j = i; j >= 0;) {
                    auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
                    if (pos == 0) {
                        continue;
                    }
                    if (indel.count(pos) == 0) {
                        indel[pos] = ALN_GAP;
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
        if (indel.empty()) {
            seq_read = seq_ref;
            return;
        }
        int k = 0;
        for (auto i = 0; i < seq_ref.size();) {
            if (indel.find(k) == indel.end()) {
                seq_read.push_back(seq_ref[i]);
                i++;
                k++;
            } else if (indel[k] == ALN_GAP) {
                i++;
                k++;
            } else {
                seq_read.push_back(indel[k]);
                k++;
            }
        }
        while (indel.find(k) != indel.end()) {
            seq_read.push_back(indel[k]);
            k++;
        }
    }

    ArtRead::ArtRead(const ArtParams& art_params, Rprob& rprob)
        : art_params_(art_params)
        , rprob_(rprob)
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