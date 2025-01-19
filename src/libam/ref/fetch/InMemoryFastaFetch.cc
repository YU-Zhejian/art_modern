#include "InMemoryFastaFetch.hh"

#include <cstddef>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "htslib/hts.h"

#include "BaseFastaFetch.hh"
#include "libam/ref/parser/fasta_parser.hh"

namespace labw::art_modern {
namespace {
    std::vector<hts_pos_t> get_seq_lengths(const std::vector<std::string>& seq)
    {
        std::vector<hts_pos_t> retv;
        retv.reserve(seq.size());
        for (auto const& s : seq) {
            retv.emplace_back(s.size());
        }
        return retv;
    }

    std::tuple<std::vector<std::string>, std::vector<std::string>> get_seq_map(const std::string& file_name)
    {
        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;

        auto file_reader = std::ifstream(file_name);
        FastaIterator fai(file_reader);
        while (true) {
            try {
                auto [id, sequence] = fai.next();
                seq_names.emplace_back(std::move(id));
                seqs.emplace_back(std::move(sequence));
            } catch (EOFException&) {
                break;
            }
        }
        file_reader.close();
        return { std::move(seq_names), std::move(seqs) };
    }
} // namespace

InMemoryFastaFetch::InMemoryFastaFetch() = default;

InMemoryFastaFetch::InMemoryFastaFetch(const std::string& file_name)
    : InMemoryFastaFetch(get_seq_map(file_name))
{
}

std::string InMemoryFastaFetch::fetch(const size_t seq_id, const hts_pos_t start, const hts_pos_t end)
{
    return seqs_[seq_id].substr(start, end - start);
}

InMemoryFastaFetch::InMemoryFastaFetch(std::vector<std::string> seq_name, std::vector<std::string> seq)
    : BaseFastaFetch(std::move(seq_name), get_seq_lengths(seq))
    , seqs_(std::move(seq))
{
}
InMemoryFastaFetch::InMemoryFastaFetch(std::tuple<std::vector<std::string>, std::vector<std::string>> seq_map)
    : BaseFastaFetch(std::move(std::get<0>(seq_map)), get_seq_lengths(std::get<1>(seq_map)))
    , seqs_(std::move(std::get<1>(seq_map)))
{
}
std::string InMemoryFastaFetch::fetch(const size_t seq_id) { return seqs_[seq_id]; }

InMemoryFastaFetch::InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept
    : InMemoryFastaFetch(std::move(other.seq_names_), std::move(other.seqs_))
{
}
InMemoryFastaFetch::InMemoryFastaFetch(
    const InMemoryFastaFetch& other, const std::ptrdiff_t from, const std::ptrdiff_t to)
    : InMemoryFastaFetch({ other.seq_names_.begin() + from, other.seq_names_.begin() + to },
          { other.seqs_.begin() + from, other.seqs_.begin() + to })
{
}

InMemoryFastaFetch::~InMemoryFastaFetch() = default;
} // namespace labw::art_modern