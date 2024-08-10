#include "seq_utils.hh"
#include "PairwiseAlignment.hh"

#include <algorithm>
#include <map>
#include <vector>

using namespace std;
namespace labw {
namespace art_modern {
    const map<char, char>& rev_comp_trans {
        { 'A', 'T' }, { 'T', 'A' }, { 'C', 'G' }, { 'G', 'C' }
    };

    std::string qual_to_str(const std::vector<int>& qual)
    {
        string retq;
        retq.resize(qual.size());
        for (size_t k = 0; k < qual.size(); k++) {
            retq[k] = (char)(qual[k] + 33);
        }
        return retq;
    }

    std::string comp(const std::string& dna)
    {
        std::string rets;
        rets.reserve(dna.length());
        for (auto i : dna) {
            rets += rev_comp_trans.at(i);
        }
        return rets;
    }
    std::string revcomp(const std::string& dna)
    {
        auto rets = comp(dna);
        reverse(rets.begin(), rets.end());
        return rets;
    }
} // namespace art_modern
} // namespace labw