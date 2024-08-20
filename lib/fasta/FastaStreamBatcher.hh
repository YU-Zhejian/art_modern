#pragma once
#include "fasta/InMemoryFastaFetch.hh"
#include "fasta/fasta_parser.hh"

#include <memory>
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

        FastaStreamBatcher(int batch_size, std::istream& stream);

    private:
        int batch_size_;
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

        InMemoryFastaStreamBatcher(int batch_size, std::shared_ptr<BaseFastaFetch> stream);

    private:
        int batch_size_;
        size_t current_index_;
        std::shared_ptr<BaseFastaFetch> stream_;
        std::mutex mutex_;
        std::vector<std::string> seq_names_;
    };

}
}