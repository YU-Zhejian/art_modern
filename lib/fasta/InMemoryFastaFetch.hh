#pragma once
#include "BaseFastaFetch.hh"
namespace labw {
namespace art_modern {

    class InMemoryFastaFetch : public BaseFastaFetch {
    public:
        explicit InMemoryFastaFetch(const std::string& file_name);
        explicit InMemoryFastaFetch(std::map<std::string, std::string, std::less<>> seq_map);
        InMemoryFastaFetch(const std::string& contig_name, const std::string& seq);
        char* cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end) override;
        std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end) override;
        ~InMemoryFastaFetch() override;

    private:
        std::map<std::string, std::string, std::less<>> seq_map_;
    };
}
}
