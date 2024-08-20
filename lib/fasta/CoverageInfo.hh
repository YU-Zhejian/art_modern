#pragma once
#include <istream>
#include <map>
#include <string>

namespace labw {
namespace art_modern {

    class CoverageInfo {
        using coverage_map = std::map<std::string, double, std::less<>>;

    public:
        explicit CoverageInfo(double static_coverage);
        CoverageInfo(coverage_map coverage_positive, coverage_map coverage_negative);
        explicit CoverageInfo(std::istream& istream);

        double coverage_positive(const std::string& contig_name) const;
        double coverage_negative(const std::string& contig_name) const;
        CoverageInfo div(int num_parts) const;

    private:
        double static_coverage_;
        coverage_map coverage_positive_;
        coverage_map coverage_negative_;
    };

} // art_modern
} // labw
