#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/ref/parser/fasta_parser.hh"

#include <htslib/hts.h>

#include <cstddef>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

    std::tuple<std::vector<std::string>, std::vector<std::string>> get_seq_map_ss(std::istream& iss)
    {
        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;
        FastaIterator fai(iss);
        while (true) {
            try {
                auto [id, sequence] = fai.next();
                seq_names.emplace_back(std::move(id));
                seqs.emplace_back(std::move(sequence));
            } catch (EOFException&) {
                break;
            }
        }
        return { std::move(seq_names), std::move(seqs) };
    }

    std::tuple<std::vector<std::string>, std::vector<std::string>> get_seq_map(const std::string& file_name)
    {
        auto file_reader = std::ifstream(file_name);
        const auto retv = get_seq_map_ss(file_reader);
        file_reader.close();
        return retv;
    }
} // namespace

InMemoryFastaFetch::InMemoryFastaFetch(const std::string& file_name)
    : InMemoryFastaFetch(get_seq_map(file_name))
{
}

InMemoryFastaFetch::InMemoryFastaFetch(std::istream& iss)
    : InMemoryFastaFetch(get_seq_map_ss(iss))
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

InMemoryFastaFetch::InMemoryFastaFetch(std::vector<std::string>&& seq_name, std::vector<std::string>&& seq)
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

InMemoryFastaFetch::InMemoryFastaFetch(
    const InMemoryFastaFetch& other, const std::ptrdiff_t from, const std::ptrdiff_t to)
    : InMemoryFastaFetch({ other.seq_names_.begin() + from, other.seq_names_.begin() + to },
          { other.seqs_.begin() + from, other.seqs_.begin() + to })
{
}

} // namespace labw::art_modern