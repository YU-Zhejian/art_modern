#pragma once
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"

#include <mutex>
#include <sstream>
#include <string>

namespace labw {
namespace art_modern {
    class FastaStreamBatcher {
    public:
        /**
         *
         * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
         * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
         */
        InMemoryFastaFetch fetch();

        FastaStreamBatcher(size_t batch_size, std::istream& stream);

    private:
        size_t batch_size_;
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

        InMemoryFastaStreamBatcher(int batch_size, BaseFastaFetch* stream);

    private:
        size_t batch_size_;
        size_t current_index_;
        BaseFastaFetch* stream_;
        std::mutex mutex_;
    };

}
}