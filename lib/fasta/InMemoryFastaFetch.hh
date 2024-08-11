//
// Created by yuzj on 24-8-11.
//

#ifndef ART_MODERN_LIB_FASTA_INMEMORYFASTAFETCH_HH
#define ART_MODERN_LIB_FASTA_INMEMORYFASTAFETCH_HH
#include "FastaFetch.hh"
namespace labw {
namespace art_modern {

    class InMemoryFastaFetch : public FastaFetch {
    public:
        explicit InMemoryFastaFetch(const std::string& file_name);
        explicit InMemoryFastaFetch(std::map<std::string, std::string, std::less<>> seq_map);
        char* cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end) override;
        std::string fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end) override;
        ~InMemoryFastaFetch() override;

    private:
        std::map<std::string, std::string, std::less<>> seq_map_;
    };
}
}

#endif // ART_MODERN_LIB_FASTA_INMEMORYFASTAFETCH_HH
