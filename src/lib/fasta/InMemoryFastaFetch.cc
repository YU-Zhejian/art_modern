#include <fstream>

#include "InMemoryFastaFetch.hh"
#include "fasta_parser.hh"

namespace labw::art_modern {

std::vector<hts_pos_t> get_seq_lengths(const std::vector<std::string>& seq)
{
    std::vector<hts_pos_t> retv;
    retv.reserve(seq.size());
    for (auto const& s : seq) {
        retv.emplace_back(s.size());
    }
    return retv;
}

InMemoryFastaFetch::InMemoryFastaFetch()
    : BaseFastaFetch() {};

std::tuple<std::vector<std::string>, std::vector<std::string>> get_seq_map(const std::string& file_name)
{
    std::vector<std::string> seq_names;
    std::vector<std::string> seqs;

    auto file_reader = std::ifstream(file_name);
    FastaIterator fai(file_reader);
    while (true) {
        try {
            auto fasta_record = fai.next();
            seq_names.emplace_back(fasta_record.id);
            seqs.emplace_back(fasta_record.sequence);
        } catch (EOFException&) {
            break;
        }
    }
    file_reader.close();
    return std::tie(seq_names, seqs);
}

InMemoryFastaFetch::InMemoryFastaFetch(const std::string& file_name)
    : InMemoryFastaFetch(get_seq_map(file_name))
{
}

InMemoryFastaFetch::InMemoryFastaFetch(const std::string& contig_name, const std::string& seq)
    : InMemoryFastaFetch(std::vector<std::string> { { contig_name } }, std::vector<std::string> { { seq } })
{
}
std::string InMemoryFastaFetch::fetch(const size_t seq_id, const hts_pos_t start, const hts_pos_t end)
{
    return seqs_[seq_id].substr(start, end - start);
}

InMemoryFastaFetch::InMemoryFastaFetch(const std::vector<std::string>& seq_name, const std::vector<std::string>& seq)
    : BaseFastaFetch(seq_name, get_seq_lengths(seq))
    , seqs_(seq)
{
}
InMemoryFastaFetch::InMemoryFastaFetch(const std::tuple<std::vector<std::string>, std::vector<std::string>>& seq_map)
    : BaseFastaFetch(std::get<0>(seq_map), get_seq_lengths(std::get<1>(seq_map)))
    , seqs_(std::get<1>(seq_map))
{
}
std::string InMemoryFastaFetch::fetch(const size_t seq_id) { return seqs_[seq_id]; }

InMemoryFastaFetch::InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept
    : InMemoryFastaFetch(other.seq_names_, other.seqs_)
{
}

InMemoryFastaFetch::~InMemoryFastaFetch() = default;
}