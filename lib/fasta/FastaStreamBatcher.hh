#pragma once
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"

#include <mutex>
#include <sstream>
#include <string>

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

    FastaStreamBatcher(const FastaStreamBatcher&) = delete;
    FastaStreamBatcher(FastaStreamBatcher&&) = delete;
    FastaStreamBatcher& operator=(const FastaStreamBatcher&) = delete;
    FastaStreamBatcher& operator=(FastaStreamBatcher&&) = delete;

private:
    std::size_t batch_size_;
    FastaIterator fasta_iterator_;
    std::mutex mutex_;
};
class InMemoryFastaStreamBatcher {
public:
    /**
     *
     * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
     * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
     */
    InMemoryFastaFetch fetch();

    InMemoryFastaStreamBatcher(std::size_t batch_size, BaseFastaFetch* stream);

    InMemoryFastaStreamBatcher(const InMemoryFastaStreamBatcher&) = delete;
    InMemoryFastaStreamBatcher(InMemoryFastaStreamBatcher&&) = delete;
    InMemoryFastaStreamBatcher& operator=(const InMemoryFastaStreamBatcher&) = delete;
    InMemoryFastaStreamBatcher& operator=(InMemoryFastaStreamBatcher&&) = delete;

private:
    std::size_t batch_size_;
    std::size_t current_index_;
    BaseFastaFetch* stream_;
    std::mutex mutex_;
};

}