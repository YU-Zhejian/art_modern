#pragma once

#include "ref/fetch/BaseFastaFetch.hh"

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include <htslib/hts.h>

namespace labw::art_modern {

/**
 * This class has move constructor only.
 */
class InMemoryFastaFetch : public BaseFastaFetch {
public:
    InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept;
    InMemoryFastaFetch(const InMemoryFastaFetch& other) = delete;
    InMemoryFastaFetch(const InMemoryFastaFetch& other, std::ptrdiff_t from, std::ptrdiff_t to);
    InMemoryFastaFetch();
    explicit InMemoryFastaFetch(const std::string& file_name);
    explicit InMemoryFastaFetch(std::tuple<std::vector<std::string>, std::vector<std::string>> seq_map);
    InMemoryFastaFetch(std::vector<std::string> seq_name, std::vector<std::string> seq);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    std::string fetch(size_t seq_id) override;
    ~InMemoryFastaFetch() override;

    InMemoryFastaFetch& operator=(const InMemoryFastaFetch&) = delete;
    InMemoryFastaFetch& operator=(InMemoryFastaFetch&&) = delete;

private:
    std::vector<std::string> seqs_;
};
} // namespace labw::art_modern
