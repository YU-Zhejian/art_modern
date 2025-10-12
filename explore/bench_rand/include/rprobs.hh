#pragma once

#include "class_utils.hh"

#include <absl/base/attributes.h>

#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <string>
#include <thread>

constexpr std::size_t N_BASES = (1 << 10);
constexpr std::size_t N_TIMES = (1 << 10);
constexpr std::size_t N_REPLICA = 200UL;
constexpr std::size_t a = 0;
constexpr std::size_t b = 1000;

ABSL_ATTRIBUTE_ALWAYS_INLINE inline static uint64_t seed()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
               .count()
        * static_cast<uint64_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));
}

static std::string formatWithCommas(std::size_t number)
{
    std::string numStr = std::to_string(number);
    int insertPosition = static_cast<int>(numStr.length()) - 3;

    while (insertPosition > 0) {
        numStr.insert(insertPosition, ",");
        insertPosition -= 3;
    }

    return numStr;
}

class DumbRandomDevice {
public:
    DELETE_COPY_MOVE(DumbRandomDevice)
    using result_type = std::uint32_t;

    DumbRandomDevice() = default;
    ~DumbRandomDevice() = default;

    SPAN_RESULT_TYPE

    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type operator()() { return 0; }
};

class CustomRandomDevice {
public:
    DELETE_COPY_MOVE(CustomRandomDevice)
    using result_type = std::uint32_t;

    explicit CustomRandomDevice(const std::string& path = "/dev/random")
        : rand_stream_(path, std::ios::binary)
        , path_(path)
    {
        if (!rand_stream_.is_open()) {
            throw std::system_error(errno, std::generic_category(), "Failed to open " + path_);
        }
    }
    ~CustomRandomDevice() { rand_stream_.close(); }

    SPAN_RESULT_TYPE

    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type operator()()
    {
        result_type random_number = 0;
        rand_stream_.read(reinterpret_cast<char*>(&random_number), sizeof(random_number));
        if (rand_stream_.fail()) {
            throw std::system_error(errno, std::generic_category(), "Failed to read from " + path_);
        }
        return random_number;
    }

private:
    std::ifstream rand_stream_;
    std::string path_;
};
