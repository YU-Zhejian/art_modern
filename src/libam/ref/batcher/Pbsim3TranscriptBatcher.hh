#pragma once

#include "libam/ds/CoverageInfo.hh"
#include "libam/ref/fetch/InMemoryFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <cstddef>
#include <istream>
#include <mutex>
#include <utility>

namespace labw::art_modern {

class Pbsim3TranscriptBatcher {
public:
    DELETE_COPY(Pbsim3TranscriptBatcher)
    DELETE_MOVE(Pbsim3TranscriptBatcher)
    ~Pbsim3TranscriptBatcher() = default;

    explicit Pbsim3TranscriptBatcher(std::size_t batch_size, std::istream& istream);
    std::pair<InMemoryFastaFetch, CoverageInfo> fetch();

private:
    std::size_t batch_size_;
    std::istream& istream_;
    std::mutex mutex_;
};

} // namespace labw::art_modern
