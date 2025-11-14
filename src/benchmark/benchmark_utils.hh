#pragma once

#include "libam_support/utils/si_utils.hh"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace labw::art_modern {
template <typename T> T geometric_mean(const std::vector<T>& data)
{
    double log_sum = 0.0;
    for (const auto& value : data) {
        log_sum += std::log(value);
    }
    return static_cast<T>(std::exp(log_sum / static_cast<double>(data.size())));
}

template <typename T> T mean(const std::vector<T>& data)
{
    return std::accumulate(data.begin(), data.end(), static_cast<T>(0)) / static_cast<T>(data.size());
}

template <typename T> T sd(const std::vector<T>& data, T mean_)
{
    T sum_squared_diff = 0;
    for (const auto& value : data) {
        T diff = value - mean_;
        sum_squared_diff += diff * diff;
    }
    return std::sqrt(sum_squared_diff / static_cast<T>(data.size() - 1));
}

template <typename T> static std::string describe(const std::vector<T>& times, bool commas = true)
{
    auto const mean_ = mean(times);
    std::ostringstream oss;
    if (commas) {
        oss << "gmean: " << std::setw(10) << format_with_commas(geometric_mean(times)) << "; mean/sd: " << std::setw(15)
            << format_with_commas(mean_) + "/" + format_with_commas(sd(times, mean_));
    } else {
        oss << "gmean: " << std::setw(10) << geometric_mean(times) << "; mean/sd: " << std::setw(15)
            << std::to_string(mean_) + "/" + std::to_string(sd(times, mean_));
    }
    return oss.str();
}

} // namespace labw::art_modern
