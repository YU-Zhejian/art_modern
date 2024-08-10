#include "ArtRead.hh"
#include "PairwiseAlignment.hh"
#include "art_modern_constants.hh"
#include "misc.hh"

#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>

using namespace std;
using namespace labw::art_modern;

namespace labw {
namespace art_modern {
    char rand_base()
    {
        static boost::random::mt19937 eng;
        boost::random::uniform_smallint<> rng(0, 3);
        switch (rng(eng)) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        default:
            return 'T';
        }
    }
} // namespace art_modern
}; // namespace labw

PairwiseAlignment ArtRead::generate_pairwise_aln() const
{
    std::string aln_read = seq_read;
    std::string aln_ref = seq_ref;
    for (auto& it : indel) {
        (it.second != ALN_GAP ? aln_ref : aln_read).insert(it.first, 1, ALN_GAP);
    }
    return { aln_read, aln_ref };
}

std::vector<int> ArtRead::generate_snv_on_qual(const std::vector<int>& qual)
{
    auto qual_mutated = qual;
    for (auto i = 0; i < qual.size(); i++) {
        if (seq_read[i] == 'N') {
            qual_mutated[i] = 1; // FIXME: This is the reason why this function cannot be const
            continue;
        }
        if (r_prob() < _art_params.err_prob[qual_mutated[i]]) {
            char achar = seq_read[i];
            while (seq_read[i] == achar) {
                achar = rand_base();
            }
            seq_read[i] = achar;
            substitution[i] = achar;
        }
    }
    return qual_mutated;
}

int ArtRead::generate_indels(int read_len, bool is_read_1)
{
    int ins_len = 0;
    int del_len = 0;
    indel.clear();
    auto per_base_del_rate = is_read_1 ? _art_params.per_base_del_rate_1
                                       : _art_params.per_base_del_rate_2;
    auto per_base_ins_rate = is_read_1 ? _art_params.per_base_ins_rate_1
                                       : _art_params.per_base_ins_rate_2;
    // deletion
    for (int i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= r_prob()) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = static_cast<int>(
                    floor((read_len - 1) * r_prob())); // invalid deletion positions: 0 or read_len-1
                if (indel.count(pos) == 0) {
                    indel[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }

    for (auto i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if ((read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= r_prob()) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = static_cast<int>(floor(r_prob() * read_len));
                if (indel.count(pos) == 0) {
                    indel[pos] = rand_base();
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len - del_len);
};

// number of deletions <= number of insertions
int ArtRead::generate_indels_2(int read_len, bool is_read_1)
{
    int ins_len = 0;
    int del_len = 0;
    indel.clear();
    auto per_base_del_rate = is_read_1 ? _art_params.per_base_del_rate_1
                                       : _art_params.per_base_del_rate_2;
    auto per_base_ins_rate = is_read_1 ? _art_params.per_base_ins_rate_1
                                       : _art_params.per_base_ins_rate_2;

    for (int i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= r_prob()) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = static_cast<int>(floor(r_prob() * read_len));
                if (indel.count(pos) == 0) {
                    indel[pos] = rand_base();
                    j--;
                }
            }
            break;
        }
    }

    // deletion
    for (int i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if ((read_len - del_len - ins_len) < (i + 1))
            continue; // ensure that enough unchanged position for mutation

        if (per_base_del_rate[i] >= r_prob()) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = static_cast<int>(
                    floor((read_len - 1) * r_prob())); // invalid deletion positions: 0 or read_len-1
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
    // I suppose this error exists no more.
    //    if (seq_read.size() != _art_params.read_len) {
    //        std::cout << map_to_str(indel) << std::endl;
    //        throw std::exception();
    //    }
}
