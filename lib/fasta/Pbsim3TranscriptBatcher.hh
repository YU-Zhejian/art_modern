#pragma once
#include "fasta/CoverageInfo.hh"
#include "fasta/InMemoryFastaFetch.hh"
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

namespace labw {
namespace art_modern {

    class Pbsim3TranscriptBatcher {
    public:
        explicit Pbsim3TranscriptBatcher(int batch_size, std::istream& istream);
        std::pair<InMemoryFastaFetch, CoverageInfo> fetch();

    private:
        size_t batch_size_;
        std::istream& istream_;
        std::mutex mutex_;
    };

} // art_modern
} // labw
