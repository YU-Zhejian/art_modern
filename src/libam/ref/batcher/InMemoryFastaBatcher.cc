#include "InMemoryFastaBatcher.hh"

#include "libam/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <cstddef>
#include <mutex>

namespace labw::art_modern {

InMemoryFastaFetch InMemoryFastaBatcher::fetch()
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

InMemoryFastaBatcher::InMemoryFastaBatcher(const std::size_t batch_size, const InMemoryFastaFetch& stream)
    : batch_size_(batch_size)
    , current_index_(0)
    , stream_(stream)
{
}
} // namespace labw::art_modern
