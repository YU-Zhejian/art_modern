#pragma once
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"

#include <mutex>
#include <sstream>
#include <string>

namespace labw {
namespace art_modern {
    class FastaStreamBatcher {
        explicit FastaStreamBatcher(int batch_size, std::istream& stream);
        /**
         *
         * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
         * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
         */
        InMemoryFastaFetch fetch();

    private:
        int batch_size_;
        FastaIterator fasta_iterator_;
        std::mutex mutex_;
    };

}
}