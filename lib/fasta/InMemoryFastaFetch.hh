#pragma once
#include "BaseFastaFetch.hh"
namespace labw {
namespace art_modern {

    class InMemoryFastaFetch : public BaseFastaFetch {
    public:
        InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept
            : InMemoryFastaFetch(other.seq_names_, other.seqs_)
        {
        }

        InMemoryFastaFetch(const InMemoryFastaFetch& other) noexcept
            : InMemoryFastaFetch(other.seq_names_, other.seqs_)
        {
        }

        InMemoryFastaFetch();
        explicit InMemoryFastaFetch(const std::string& file_name);
        explicit InMemoryFastaFetch(const std::unordered_map<std::string, std::string>& seq_map);
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
}
