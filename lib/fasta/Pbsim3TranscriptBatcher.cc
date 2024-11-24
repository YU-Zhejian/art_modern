#include "Pbsim3TranscriptBatcher.hh"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <utility>
#include <vector>

namespace labw {
namespace art_modern {
    Pbsim3TranscriptBatcher::Pbsim3TranscriptBatcher(std::size_t batch_size, std::istream& istream)
        : batch_size_(batch_size)
        , istream_(istream)
    {
    }
    std::pair<InMemoryFastaFetch, CoverageInfo> Pbsim3TranscriptBatcher::fetch()
    {
        std::scoped_lock lock(mutex_);
        std::unordered_map<std::string, double> coverage_positive;
        std::unordered_map<std::string, double> coverage_negative;
        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;
        if (batch_size_ != std::numeric_limits<int>::max()) {
            coverage_positive.reserve(batch_size_);
            coverage_negative.reserve(batch_size_);

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
            boost::algorithm::split(tokens, line, boost::is_any_of("\t"));
            if (tokens.size() == 4) {
                coverage_positive[tokens.at(0)] = std::stod(tokens.at(1));
                coverage_negative[tokens.at(0)] = std::stod(tokens.at(2));
                seq_names.emplace_back(tokens.at(0));
                seqs.emplace_back(tokens.at(3));
            } else {
                throw std::invalid_argument("Cannot parse PBSIM3 transcript " + line);
            }
        }
        return { InMemoryFastaFetch(seq_names, seqs), CoverageInfo(coverage_positive, coverage_negative) };
    }
} // art_modern
} // labw