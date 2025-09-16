#include "InMemoryFastaBatcher.hh"

#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mutex>

namespace labw::art_modern {
std::shared_ptr<InMemoryFastaFetch> InMemoryFastaBatcher::fetch()
{
    const std::scoped_lock lock(mutex_);
    const auto from = current_index_;
    const auto to = std::min(current_index_ + batch_size_, stream_->num_seqs());
    auto fetch
        = std::make_shared<InMemoryFastaFetch>(*stream_, static_cast<ptrdiff_t>(from), static_cast<ptrdiff_t>(to));
    if (!fetch->empty()) {
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch->seq_name(0) << " to "
                                << fetch->seq_name(fetch->num_seqs() - 1) << " (" << fetch->num_seqs() << "） created";
    }
    current_index_ = to;
    return fetch;
}

InMemoryFastaBatcher::InMemoryFastaBatcher(
    const std::size_t batch_size, const std::shared_ptr<InMemoryFastaFetch>& stream)
    : batch_size_(batch_size)
    , stream_(stream)
{
}
} // namespace labw::art_modern
