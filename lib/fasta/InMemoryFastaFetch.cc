#include <fstream>

#include "InMemoryFastaFetch.hh"
#include "fasta_parser.hh"

namespace labw {
namespace art_modern {
    InMemoryFastaFetch::InMemoryFastaFetch(std::map<std::string, std::string, std::less<>> seq_map)
        : seq_map_(std::move(seq_map))
    {
        for (auto const& pair : seq_map_) {
            seq_lengths_.emplace(pair.first, pair.second.size());
        }
    }
    std::string InMemoryFastaFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        return seq_map_.at(seq_name).substr(start, end - start);
    }
    char* InMemoryFastaFetch::cfetch(const char* seq_name, hts_pos_t start, hts_pos_t end)
    {
        auto fetch_str = fetch(seq_name, start, end);
        auto rets = (char*)calloc(fetch_str.size() + 1, sizeof(char));
        strncpy(rets, fetch_str.c_str(), fetch_str.size());
        return rets;
    }
    InMemoryFastaFetch::InMemoryFastaFetch(const std::string& file_name)
    {
        auto file_reader = std::ifstream(file_name);
        FastaIterator fai(file_reader);
        while (true) {
            try {
                auto fasta_record = fai.next();
                seq_map_.emplace(fasta_record.id, fasta_record.sequence);
                seq_lengths_.emplace(fasta_record.id, fasta_record.sequence.size());
            } catch (EOFException&) {
                break;
            }
        }
        file_reader.close();
    }

    InMemoryFastaFetch::~InMemoryFastaFetch() = default;
}
}