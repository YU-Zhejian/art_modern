#include <boost/log/trivial.hpp>

#include "BaseFastaFetch.hh"
#include "NotImplementedException.hh"

namespace labw {
namespace art_modern {
    std::string BaseFastaFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        throw NotImplementedException();
    }
    char* BaseFastaFetch::cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        throw NotImplementedException();
    }
    void BaseFastaFetch::update_sam_header(sam_hdr_t* header) const
    {
        for (const auto& pair : seq_lengths_) {
            auto seq_len_str = std::to_string(pair.second);
            sam_hdr_add_line(header, "SQ", "SN", pair.first.c_str(), "LN", seq_len_str.c_str(), NULL);
        }
    }
    hts_pos_t BaseFastaFetch::seq_len(const std::string& seq_name) const
    {
        return seq_lengths_.at(seq_name);
    }
    size_t BaseFastaFetch::num_seqs() const
    {
        return seq_lengths_.size();
    }
    BaseFastaFetch::~BaseFastaFetch() = default;
}
}
