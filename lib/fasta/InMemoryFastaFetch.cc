#include <fstream>

#include "InMemoryFastaFetch.hh"
#include "MapUtils.hh"
#include "fasta_parser.hh"

namespace labw {
namespace art_modern {

    std::vector<hts_pos_t> get_seq_lengths(const std::vector<std::string>& seq)
    {
        std::vector<hts_pos_t> retv;
        retv.reserve(seq.size());
        for (auto const& s : seq) {
            retv.emplace_back(s.size());
        }
        return retv;
    }

    InMemoryFastaFetch::InMemoryFastaFetch(const std::unordered_map<std::string, std::string>& seq_map)
        : InMemoryFastaFetch(convert_map_to_k_v_list(seq_map))
    {
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
    std::string InMemoryFastaFetch::fetch(size_t seq_id, hts_pos_t start, hts_pos_t end)
    {
        return seqs_[seq_id].substr(start, end - start);
    }

    InMemoryFastaFetch::InMemoryFastaFetch(
        const std::vector<std::string>& seq_name, const std::vector<std::string>& seq)
        : BaseFastaFetch(seq_name, get_seq_lengths(seq))
        , seqs_(seq)
    {
    }
    InMemoryFastaFetch::InMemoryFastaFetch(
        const std::tuple<std::vector<std::string>, std::vector<std::string>>& seq_map)
        : BaseFastaFetch(std::get<0>(seq_map), get_seq_lengths(std::get<1>(seq_map)))
        , seqs_(std::get<1>(seq_map))
    {
    }
    std::string InMemoryFastaFetch::fetch(size_t seq_id) { return seqs_[seq_id]; }

    InMemoryFastaFetch::~InMemoryFastaFetch() = default;
}
}