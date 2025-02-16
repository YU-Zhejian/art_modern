#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <istream>
#include <string>
#include <tuple>
#include <unordered_map>

namespace labw::art_modern {

class CoverageInfo {

public:
    using coverage_map = std::unordered_map<std::string, double>;

    explicit CoverageInfo(double static_coverage);
    explicit CoverageInfo(double static_coverage_positive, double static_coverage_negative);
    CoverageInfo(coverage_map&& coverage_positive, coverage_map&& coverage_negative);
    explicit CoverageInfo(std::istream& istream);
    explicit CoverageInfo(const std::tuple<coverage_map, coverage_map>& coverage);

    ~CoverageInfo() = default;

    DEFAULT_COPY(CoverageInfo)
    DEFAULT_MOVE(CoverageInfo)

    [[nodiscard]] double coverage_positive(const std::string& contig_name) const;
    [[nodiscard]] double coverage_negative(const std::string& contig_name) const;
    [[nodiscard]] CoverageInfo div(int num_parts) const;

private:
    double static_coverage_positive_;
    double static_coverage_negative_;
    coverage_map coverage_positive_;
    coverage_map coverage_negative_;
};

} // namespace labw::art_modern
