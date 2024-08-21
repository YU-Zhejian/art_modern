#pragma once
#include "BaseFastaFetch.hh"
namespace labw {
namespace art_modern {

    class InMemoryFastaFetch : public BaseFastaFetch {
    public:
        InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept: InMemoryFastaFetch(other.seq_map_){
        }

        InMemoryFastaFetch(const InMemoryFastaFetch & other) noexcept:
            InMemoryFastaFetch(other.seq_map_){
        }

        InMemoryFastaFetch();
        explicit InMemoryFastaFetch(const std::string& file_name);
        explicit InMemoryFastaFetch(std::map<std::string, std::string, std::less<>> seq_map);
        InMemoryFastaFetch(const std::string& contig_name, const std::string& seq);
        std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end) override;
        ~InMemoryFastaFetch() override;

    private:
        const std::map<std::string, std::string, std::less<>> seq_map_;
    };
}
}
