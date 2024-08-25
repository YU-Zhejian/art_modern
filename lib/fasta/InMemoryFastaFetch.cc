#include <fstream>

#include "InMemoryFastaFetch.hh"
#include "fasta_parser.hh"

namespace labw {
namespace art_modern {

    std::unordered_map<std::string, hts_pos_t> get_seq_lengths(
        const std::unordered_map<std::string, std::string>& seq_map)
    {
        std::unordered_map<std::string, hts_pos_t> seq_lengths;

        for (auto const& pair : seq_map) {
            seq_lengths.emplace(pair.first, pair.second.size());
        }
        return seq_lengths;
    }

    InMemoryFastaFetch::InMemoryFastaFetch(std::unordered_map<std::string, std::string> seq_map)
        : BaseFastaFetch(get_seq_lengths(seq_map))
        , seq_map_(std::move(seq_map))
    {
    }

    std::string InMemoryFastaFetch::fetch(const std::string& seq_name, hts_pos_t start, hts_pos_t end)
    {
        return seq_map_.at(seq_name).substr(start, end - start);
    }

    InMemoryFastaFetch::InMemoryFastaFetch()
        : BaseFastaFetch(std::unordered_map<std::string, hts_pos_t>()) {};

    std::unordered_map<std::string, std::string> get_seq_map(const std::string& file_name)
    {
        std::unordered_map<std::string, std::string> seq_map;
        auto file_reader = std::ifstream(file_name);
        FastaIterator fai(file_reader);
        while (true) {
            try {
                auto fasta_record = fai.next();
                seq_map.emplace(fasta_record.id, fasta_record.sequence);
            } catch (EOFException&) {
                break;
            }
        }
        file_reader.close();
        return seq_map;
    }

    InMemoryFastaFetch::InMemoryFastaFetch(const std::string& file_name)
        : InMemoryFastaFetch(get_seq_map(file_name))
    {
    }

    InMemoryFastaFetch::InMemoryFastaFetch(const std::string& contig_name, const std::string& seq)
        : InMemoryFastaFetch(std::unordered_map<std::string, std::string> { { contig_name, seq } })
    {
    }

    InMemoryFastaFetch::~InMemoryFastaFetch() = default;
}
}