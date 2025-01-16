#pragma once

#include "libam/ref/fetch/InMemoryFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <cstddef>
#include <mutex>

namespace labw::art_modern {

class InMemoryFastaBatcher {
public:
    /**
     *
     * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
     * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
     */
    InMemoryFastaFetch fetch();

    InMemoryFastaBatcher(std::size_t batch_size, const InMemoryFastaFetch& stream);
    ~InMemoryFastaBatcher() = default;

    DELETE_COPY(InMemoryFastaBatcher)
    DELETE_MOVE(InMemoryFastaBatcher)

private:
    std::size_t batch_size_;
    std::size_t current_index_;
    const InMemoryFastaFetch& stream_;
    std::mutex mutex_;
};

} // namespace labw::art_modern
