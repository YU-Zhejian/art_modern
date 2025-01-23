#include "libam/ref/fetch/BaseFastaFetch.hh"

#include "art_modern_config.h"
#include "libam/CExceptionsProxy.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace labw::art_modern {

void BaseFastaFetch::update_sam_header(sam_hdr_t* header) const
{
    for (size_t i = 0; i < seq_lengths_.size(); ++i) {
        auto seq_len_str = std::to_string(seq_lengths_[i]);
        CExceptionsProxy::assert_numeric(
            sam_hdr_add_line(header, "SQ", "SN", seq_names_[i].c_str(), "LN", seq_len_str.c_str(), NULL),
            USED_HTSLIB_NAME, "Failed to add SQ header line for contig '" + seq_names_[i] + "'", false,
            CExceptionsProxy::EXPECTATION::ZERO);
    }
}

hts_pos_t BaseFastaFetch::seq_len(const std::size_t seq_id) const { return seq_lengths_[seq_id]; }

std::string BaseFastaFetch::seq_name(const std::size_t seq_id) const { return seq_names_[seq_id]; }

size_t BaseFastaFetch::num_seqs() const { return seq_lengths_.size(); }

BaseFastaFetch::BaseFastaFetch(std::vector<std::string>&& seq_names, std::vector<hts_pos_t>&& seq_lengths)
    : seq_names_(std::move(seq_names))
    , seq_lengths_(std::move(seq_lengths))
{
}
BaseFastaFetch::BaseFastaFetch(const std::vector<std::string>& seq_names, const std::vector<hts_pos_t>& seq_lengths)
    : seq_names_(seq_names)
    , seq_lengths_(seq_lengths)
{
}
BaseFastaFetch::BaseFastaFetch(const std::tuple<std::vector<std::string>, std::vector<hts_pos_t>>& seq_names_lengths)
    : seq_names_(std::get<0>(seq_names_lengths))
    , seq_lengths_(std::get<1>(seq_names_lengths))
{
}
bool BaseFastaFetch::empty() const { return this->seq_names_.empty(); }
std::string BaseFastaFetch::fetch(const std::size_t seq_id) { return fetch(seq_id, 0, seq_lengths_[seq_id]); }
} // namespace labw::art_modern
