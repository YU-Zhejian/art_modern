#include <boost/log/trivial.hpp>

#include "BaseFastaFetch.hh"
#include "CExceptionsProxy.hh"
#include "ExceptionUtils.hh"
#include "art_modern_constants.hh"

namespace labw {
namespace art_modern {

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
    hts_pos_t BaseFastaFetch::seq_len(const std::string& seq_name) const
    {
        try {
            return seq_lengths_.at(seq_name);
        } catch (std::out_of_range&) {
            BOOST_LOG_TRIVIAL(error) << "Invalid sequence name " << seq_name << ".";
            throw std::invalid_argument("Invalid sequence name.");
        }
    }

    size_t BaseFastaFetch::num_seqs() const { return seq_lengths_.size(); }

    const std::vector<std::string>& BaseFastaFetch::seq_names() const { return seq_names_; }

    std::vector<std::string> get_seq_names(const std::unordered_map<std::string, hts_pos_t>& seq_lengths)
    {
        std::vector<std::string> retv;
        retv.reserve(seq_lengths.size());
        for (const auto& pair : seq_lengths) {
            retv.emplace_back(pair.first);
        }
        return retv;
    }

    BaseFastaFetch::BaseFastaFetch(std::unordered_map<std::string, hts_pos_t> seq_lengths)
        : seq_lengths_(std::move(seq_lengths))
        , seq_names_(get_seq_names(seq_lengths_))
    {
    }

    BaseFastaFetch::~BaseFastaFetch() = default;
}
}
