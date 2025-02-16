#pragma once
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/ref/parser/fasta_parser.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <istream>
#include <mutex>

namespace labw::art_modern {
class FastaStreamBatcher {
public:
    /**
     *
     * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
     * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
     */
    InMemoryFastaFetch fetch();

    FastaStreamBatcher(std::size_t batch_size, std::istream& stream);

    DELETE_COPY(FastaStreamBatcher)
    DELETE_MOVE(FastaStreamBatcher)
    ~FastaStreamBatcher() = default;

private:
    std::size_t batch_size_;
    FastaIterator fasta_iterator_;
    std::mutex mutex_;
};

} // namespace labw::art_modern