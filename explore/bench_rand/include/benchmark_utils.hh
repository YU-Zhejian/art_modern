#pragma once

#include "rprobs.hh"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

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

static std::string describe(const std::vector<std::size_t>& times)
{
    auto const mean_ = mean(times);
    std::ostringstream oss;
    oss << "gmean: " << std::setw(10) << formatWithCommas(geometric_mean(times)) << "; mean/sd: " << std::setw(15)
        << formatWithCommas(mean_) + "/" + formatWithCommas(sd(times, mean_));
    return oss.str();
}
