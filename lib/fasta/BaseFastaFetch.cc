#include <boost/log/trivial.hpp>

#include "BaseFastaFetch.hh"
#include "CExceptionsProxy.hh"
#include "ExceptionUtils.hh"
#include "NotImplementedException.hh"
#include "art_modern_constants.hh"

namespace labw {
namespace art_modern {
    std::string BaseFastaFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseFastaFetch::fetch?";
        throw_with_trace(NotImplementedException());
    }
    char* BaseFastaFetch::cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseFastaFetch::cfetch?";
        throw_with_trace(NotImplementedException());
    }
    void BaseFastaFetch::update_sam_header(sam_hdr_t* header) const
    {
        for (const auto& pair : seq_lengths_) {
            auto seq_len_str = std::to_string(pair.second);
            CExceptionsProxy::requires_numeric(
                sam_hdr_add_line(header, "SQ", "SN", pair.first.c_str(), "LN", seq_len_str.c_str(), NULL),
                USED_HTSLIB_NAME, "Failed to add SQ header line for contig '" + pair.first + "'", false,
                CExceptionsProxy::EXPECTATION::ZERO);
        }
    }
    hts_pos_t BaseFastaFetch::seq_len(const std::string& seq_name) const { return seq_lengths_.at(seq_name); }
    size_t BaseFastaFetch::num_seqs() const { return seq_lengths_.size(); }

    std::vector<std::string> BaseFastaFetch::seq_names() const { return seq_names_; }

    BaseFastaFetch::~BaseFastaFetch() = default;
}
}
