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
        std::lock_guard<std::mutex> lock(mutex_);
        std::unordered_map<std::string, double> coverage_positive;
        std::unordered_map<std::string, double> coverage_negative;
        coverage_positive.reserve(batch_size_);
        coverage_negative.reserve(batch_size_);

        std::vector<std::string> seq_names;
        std::vector<std::string> seqs;
        seq_names.reserve(batch_size_);
        seqs.reserve(batch_size_);

        std::string buff;
        std::vector<std::string> tokens;
        while (seq_names.size() < batch_size_ && !istream_.eof()) {
            std::getline(istream_, buff);
            if (buff.empty() || buff.at(0) == '#') {
                continue;
            }
            boost::algorithm::split(tokens, buff, boost::is_any_of("\t"));
            if (tokens.size() == 4) {
                coverage_positive[tokens.at(0)] = std::stod(tokens.at(1));
                coverage_negative[tokens.at(0)] = std::stod(tokens.at(2));
                seq_names.emplace_back(tokens.at(0));
                seqs.emplace_back(tokens.at(3));
            } else {
                throw std::invalid_argument("Cannot parse PBSIM3 transcript " + buff);
            }
        }
        return std::make_pair(InMemoryFastaFetch(seq_names, seqs), CoverageInfo(coverage_positive, coverage_negative));
    }
} // art_modern
} // labw