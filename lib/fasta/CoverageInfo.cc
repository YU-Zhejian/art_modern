#include "CoverageInfo.hh"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <sstream>
#include <utility>
#include <vector>

namespace labw {
namespace art_modern {
    CoverageInfo::CoverageInfo(const double static_coverage)
        : static_coverage_positive_(static_coverage / 2)
        , static_coverage_negative_(static_coverage / 2)
    {
    }

    CoverageInfo::CoverageInfo(
        CoverageInfo::coverage_map coverage_positive, CoverageInfo::coverage_map coverage_negative)
        : static_coverage_positive_(0)
        , static_coverage_negative_(0)
        , coverage_positive_(std::move(coverage_positive))
        , coverage_negative_(std::move(coverage_negative))
    {
    }

    double CoverageInfo::coverage_positive(const std::string& contig_name) const
    {
        if (coverage_positive_.empty()) {
            return static_coverage_positive_;
        } else {
            if (coverage_positive_.find(contig_name) == coverage_positive_.end()) {
                return 0;
            }
            return coverage_positive_.at(contig_name);
        }
    }

    double CoverageInfo::coverage_negative(const std::string& contig_name) const
    {
        if (coverage_negative_.empty()) {
            return static_coverage_negative_;
        } else {
            if (coverage_negative_.find(contig_name) == coverage_negative_.end()) {
                return 0;
            }
            return coverage_negative_.at(contig_name);
        }
    }

    std::tuple<CoverageInfo::coverage_map, CoverageInfo::coverage_map> read(std::istream& istream)
    {
        CoverageInfo::coverage_map coverage_positive;
        CoverageInfo::coverage_map coverage_negative;
        std::string buff;
        std::vector<std::string> tokens;
        while (!istream.eof()) {
            std::getline(istream, buff);
            if (buff.empty() || buff.front() == '#') {
                continue;
            }
            boost::algorithm::split(tokens, buff, boost::is_any_of("\t"));
            if (tokens.size() == 3) {
                coverage_positive[tokens.at(0)] = std::stod(tokens.at(1));
                coverage_negative[tokens.at(0)] = std::stod(tokens.at(2));
            } else if (tokens.size() == 2) {
                coverage_positive[tokens.at(0)] = std::stod(tokens.at(1)) / 2;
                coverage_negative[tokens.at(0)] = std::stod(tokens.at(1)) / 2;
            }
        }
        return { coverage_positive, coverage_negative };
    }

    CoverageInfo::CoverageInfo(std::istream& istream)
        : CoverageInfo(read(istream))
    {
    }

    CoverageInfo CoverageInfo::div(const int num_parts) const
    {
        if (coverage_positive_.empty()) {
            return CoverageInfo(static_coverage_positive_ / num_parts, static_coverage_negative_ / num_parts);
        } else {
            coverage_map coverage_positive_new;
            for (auto const& pair : coverage_positive_) {
                coverage_positive_new[pair.first] = pair.second / num_parts;
            }
            coverage_map coverage_negative_new;
            for (auto const& pair : coverage_positive_) {
                coverage_negative_new[pair.first] = pair.second / num_parts;
            }
            return { coverage_positive_new, coverage_negative_new };
        }
    }

    CoverageInfo::CoverageInfo(const double static_coverage_positive, const double static_coverage_negative)
        : static_coverage_positive_(static_coverage_positive)
        , static_coverage_negative_(static_coverage_negative)
        , coverage_positive_()
        , coverage_negative_()
    {
    }
    CoverageInfo::CoverageInfo(const std::tuple<coverage_map, coverage_map>& coverage)
        : static_coverage_positive_(0)
        , static_coverage_negative_(0)
        , coverage_positive_(std::get<0>(coverage))
        , coverage_negative_(std::get<1>(coverage))
    {
    }

} // art_modern
} // labw