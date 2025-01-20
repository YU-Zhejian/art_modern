#pragma once

#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <htslib/hts.h>

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace labw::art_modern {

/**
 * This class has move constructor only.
 */
class InMemoryFastaFetch : public BaseFastaFetch {
public:
    InMemoryFastaFetch(const InMemoryFastaFetch& other)
        : BaseFastaFetch(other.seq_names_, other.seq_lengths_)
        , seqs_(other.seqs_)
    {
    }
    InMemoryFastaFetch& operator=(const InMemoryFastaFetch&) = delete;
    DELETE_MOVE(InMemoryFastaFetch)
    InMemoryFastaFetch() = default;

    InMemoryFastaFetch(const InMemoryFastaFetch& other, std::ptrdiff_t from, std::ptrdiff_t to);
    explicit InMemoryFastaFetch(const std::string& file_name);
    explicit InMemoryFastaFetch(std::tuple<std::vector<std::string>, std::vector<std::string>> seq_map);
    InMemoryFastaFetch(std::vector<std::string>&& seq_name, std::vector<std::string>&& seq);
    InMemoryFastaFetch(const std::vector<std::string>& seq_name, const std::vector<std::string>& seq);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    std::string fetch(size_t seq_id) override;
    ~InMemoryFastaFetch() override = default;

private:
    const std::vector<std::string> seqs_;
};
} // namespace labw::art_modern
