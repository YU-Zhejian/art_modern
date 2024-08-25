#pragma once
#include <istream>
#include <string>
#include <unordered_map>

namespace labw {
namespace art_modern {

    class CoverageInfo {
        using coverage_map = std::unordered_map<std::string, double>;

    public:
        explicit CoverageInfo(double static_coverage);
        explicit CoverageInfo(double static_coverage_positive, double static_coverage_negative);
        CoverageInfo(coverage_map coverage_positive, coverage_map coverage_negative);
        explicit CoverageInfo(std::istream& istream);

        double coverage_positive(const std::string& contig_name) const;
        double coverage_negative(const std::string& contig_name) const;
        CoverageInfo div(int num_parts) const;

    private:
        const double static_coverage_positive_;
        const double static_coverage_negative_;
        const coverage_map coverage_positive_;
        const coverage_map coverage_negative_;
    };

} // art_modern
} // labw
