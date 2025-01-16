#pragma once
#include "libam/ds/CoverageInfo.hh"
#include "libam/ref/fetch/InMemoryFastaFetch.hh"
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

namespace labw::art_modern {

class Pbsim3TranscriptBatcher {
public:
    explicit Pbsim3TranscriptBatcher(std::size_t batch_size, std::istream& istream);
    std::pair<InMemoryFastaFetch, CoverageInfo> fetch();

    Pbsim3TranscriptBatcher(const Pbsim3TranscriptBatcher&) = delete;
    Pbsim3TranscriptBatcher(Pbsim3TranscriptBatcher&&) = delete;
    Pbsim3TranscriptBatcher& operator=(const Pbsim3TranscriptBatcher&) = delete;
    Pbsim3TranscriptBatcher& operator=(Pbsim3TranscriptBatcher&&) = delete;

private:
    std::size_t batch_size_;
    std::istream& istream_;
    std::mutex mutex_;
};

} // art_modern
// labw
