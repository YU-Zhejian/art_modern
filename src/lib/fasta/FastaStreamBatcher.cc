#include "FastaStreamBatcher.hh"
#include <boost/log/trivial.hpp>

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

InMemoryFastaFetch InMemoryFastaStreamBatcher::fetch()
{
    std::scoped_lock lock(mutex_);
    const auto from = current_index_;
    const auto to = std::min(current_index_ + batch_size_, stream_.num_seqs());
    InMemoryFastaFetch fetch = { stream_, static_cast<ptrdiff_t>(from), static_cast<ptrdiff_t>(to) };
    if (!fetch.empty()) {
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch.seq_name(0) << " to "
                                << fetch.seq_name(fetch.num_seqs() - 1) << " (" << fetch.num_seqs() << "） created";
    }
    current_index_ = to;
    return fetch;
}

InMemoryFastaStreamBatcher::InMemoryFastaStreamBatcher(const std::size_t batch_size, const InMemoryFastaFetch& stream)
    : batch_size_(batch_size)
    , current_index_(0)
    , stream_(stream)
{
}
}