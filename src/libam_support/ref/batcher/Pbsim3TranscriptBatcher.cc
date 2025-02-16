#include "libam_support/ref/batcher/Pbsim3TranscriptBatcher.hh"

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
Pbsim3TranscriptBatcher::Pbsim3TranscriptBatcher(const std::size_t batch_size, std::istream& istream)
    : batch_size_(batch_size)
    , istream_(istream)
{
}
std::pair<std::shared_ptr<InMemoryFastaFetch>, std::shared_ptr<CoverageInfo>> Pbsim3TranscriptBatcher::fetch()
{
    const std::scoped_lock lock(mutex_);
    CoverageInfo::coverage_map coverage_positive;
    CoverageInfo::coverage_map coverage_negative;
    std::vector<std::string> seq_names;
    std::vector<std::string> seqs;
    if (batch_size_ != std::numeric_limits<int>::max()) {
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);
    }

    std::string line;
    std::vector<std::string> tokens;
    while (seq_names.size() < batch_size_ && !istream_.eof()) {
        tokens.clear();
        std::getline(istream_, line);
        if (line.empty() || line.at(0) == '#') {
            continue;
        }
        split(tokens, line, boost::is_any_of("\t"));
        if (tokens.size() == 4) {
            coverage_positive.emplace(tokens.at(0), std::stod(tokens.at(1)));
            coverage_negative.emplace(tokens.at(0), std::stod(tokens.at(2)));
            seq_names.emplace_back(tokens.at(0));
            seqs.emplace_back(tokens.at(3));
        } else {
            throw std::invalid_argument("Cannot parse PBSIM3 transcript " + line);
        }
    }
    return { std::make_shared<InMemoryFastaFetch>(std::move(seq_names), std::move(seqs)),
        std::make_shared<CoverageInfo>(std::move(coverage_positive), std::move(coverage_negative)) };
}
} // namespace labw::art_modern
