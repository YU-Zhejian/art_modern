#include "ref/batcher/FastaStreamBatcher.hh"

#include "ref/fetch/InMemoryFastaFetch.hh"
#include "ref/parser/fasta_parser.hh"

#include <boost/log/trivial.hpp>

#include <cstddef>
#include <istream>
#include <limits>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
FastaStreamBatcher::FastaStreamBatcher(const std::size_t batch_size, std::istream& stream)
    : batch_size_(batch_size)
    , fasta_iterator_(stream)
{
}
InMemoryFastaFetch FastaStreamBatcher::fetch()
{
    std::scoped_lock lock(mutex_);
    std::vector<std::string> seq_names;
    std::vector<std::string> seqs;
    if (batch_size_ != std::numeric_limits<int>::max()) {
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);
    }

    std::string fetch_s;
    std::string fetch_e;
    while (seq_names.size() < batch_size_) {
        try {
            auto [id, sequence] = fasta_iterator_.next();
            if (fetch_s.empty()) {
                fetch_s = id;
            }
            fetch_e = id;
            seq_names.emplace_back(std::move(id));
            seqs.emplace_back(std::move(sequence));
        } catch (EOFException&) {
            break;
        }
    }
    BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch_s << " to " << fetch_e << " (" << seq_names.size()
                            << "） created";
    return { std::move(seq_names), std::move(seqs) };
}

} // namespace labw::art_modern