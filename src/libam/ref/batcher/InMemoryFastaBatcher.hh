#pragma once

#include "libam/ref/fetch/InMemoryFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <cstddef>
#include <memory>
#include <mutex>

namespace labw::art_modern {

class InMemoryFastaBatcher {
public:
    InMemoryFastaBatcher(std::size_t batch_size, const std::shared_ptr<InMemoryFastaFetch>& stream);
    ~InMemoryFastaBatcher() = default;

    DELETE_COPY(InMemoryFastaBatcher)
    DELETE_MOVE(InMemoryFastaBatcher)
    /*!
     *
     * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
     * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
     */
    std::shared_ptr<InMemoryFastaFetch> fetch();

private:
    std::size_t batch_size_;
    std::size_t current_index_;
    const std::shared_ptr<InMemoryFastaFetch>& stream_;
    std::mutex mutex_;
};

} // namespace labw::art_modern
