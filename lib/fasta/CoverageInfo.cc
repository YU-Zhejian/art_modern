#include "CoverageInfo.hh"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <utility>
#include <vector>

namespace labw {
namespace art_modern {
    CoverageInfo::CoverageInfo(const double static_coverage)
        : static_coverage_(static_coverage)
    {
    }
    CoverageInfo::CoverageInfo(
        CoverageInfo::coverage_map coverage_positive, CoverageInfo::coverage_map coverage_negative)
        : static_coverage_(0)
        , coverage_positive_(std::move(coverage_positive))
        , coverage_negative_(std::move(coverage_negative))
    {
    }
    double CoverageInfo::coverage_positive(const std::string& contig_name) const
    {
        if (coverage_positive_.empty()) {
            return static_coverage_ / 2;
        } else {
            return coverage_positive_.at(contig_name);
        }
    }
    double CoverageInfo::coverage_negative(const std::string& contig_name) const
    {
        if (coverage_negative_.empty()) {
            return static_coverage_ / 2;
        } else {
            return coverage_negative_.at(contig_name);
        }
    }
    CoverageInfo::CoverageInfo(std::istream& istream)
        : static_coverage_(0)
    {
        std::string buff;
        std::vector<std::string> tokens;
        while (!istream.eof()) {
            std::getline(istream, buff);
            if (buff.empty() || buff.at(0) == '#') {
                continue;
            }
            boost::algorithm::split(tokens, buff, boost::is_any_of("\t"));
            if (tokens.size() == 3) {
                coverage_positive_[tokens.at(0)] = std::stod(tokens.at(1));
                coverage_negative_[tokens.at(0)] = std::stod(tokens.at(2));
            } else if (tokens.size() == 2) {
                coverage_positive_[tokens.at(0)] = std::stod(tokens.at(1)) / 2;
                coverage_negative_[tokens.at(0)] = std::stod(tokens.at(1)) / 2;
            }
        }
    }
    CoverageInfo CoverageInfo::operator/=(int num_parts)
    {
        if (coverage_positive_.empty()) {
            return CoverageInfo(static_coverage_ / num_parts);
        } else {
            coverage_map coverage_positive_new;
            for (auto& pair : coverage_positive_) {
                coverage_positive_new[pair.first] = pair.second / num_parts;
            }
            coverage_map coverage_negative_new;
            for (auto& pair : coverage_positive_) {
                coverage_negative_new[pair.first] = pair.second / num_parts;
            }
            return { coverage_positive_new, coverage_negative_new };
        }
    }
} // art_modern
} // labw