#include <boost/log/trivial.hpp>

#include "FastaFetch.hh"

namespace labw {
namespace art_modern {
    std::string FastaFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        // No-op implementation
        return "";
    }
    char* FastaFetch::cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        // No-op implementation
        return nullptr;
    }
    void FastaFetch::update_sam_header(sam_hdr_t* header)
    {
        for (const auto& pair : seq_lengths_) {
            auto seq_len_str = std::to_string(pair.second);
            sam_hdr_add_line(header, "SQ", "SN", pair.first.c_str(), "LN", seq_len_str.c_str());
        }
    }
    hts_pos_t FastaFetch::seq_len(const std::string& seq_name)
    {
        return seq_lengths_.at(seq_name);
    }
    size_t FastaFetch::num_seqs()
    {
        return seq_lengths_.size();
    }
    FastaFetch::~FastaFetch() = default;
}
}
