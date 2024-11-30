#pragma once
#include "BaseFastaFetch.hh"
namespace labw::art_modern {

/**
 * This class have move constructor only.
 */
class InMemoryFastaFetch : public BaseFastaFetch {
public:
    InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept;
    InMemoryFastaFetch(const InMemoryFastaFetch& other) = delete;
    InMemoryFastaFetch();
    explicit InMemoryFastaFetch(const std::string& file_name);
    explicit InMemoryFastaFetch(const std::tuple<std::vector<std::string>, std::vector<std::string>>& seq_map);
    InMemoryFastaFetch(const std::vector<std::string>& seq_name, const std::vector<std::string>& seq);
    InMemoryFastaFetch(const std::string& contig_name, const std::string& seq);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    std::string fetch(size_t seq_id) override;
    ~InMemoryFastaFetch() override;

    InMemoryFastaFetch& operator=(const InMemoryFastaFetch&) = delete;
    InMemoryFastaFetch& operator=(InMemoryFastaFetch&&) = delete;

private:
    const std::vector<std::string> seqs_;
};
}
