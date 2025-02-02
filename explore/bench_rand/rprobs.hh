#pragma once

#include <climits>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <string>

constexpr std::size_t N_BASES = (1 << 10);
constexpr std::size_t N_TIMES = 20UL * (1 << 10);
constexpr std::size_t a = 0;
constexpr std::size_t b = 1000;

static uint64_t seed()
{
    return 0;
    //    return
    //    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
    //               .count()
    //        * static_cast<uint64_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));
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
    using result_type = unsigned char;

    DumbRandomDevice() = default;
    ~DumbRandomDevice() = default;

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return CHAR_MAX; }

    result_type operator()() { return 0; }
};

class CustomRandomDevice {
public:
    using result_type = unsigned int;

    CustomRandomDevice()
        : rand_stream("/dev/random", std::ios::binary)
    {
        if (!rand_stream.is_open()) {
            throw std::system_error(errno, std::generic_category(), "Failed to open /dev/random");
        }
    }
    ~CustomRandomDevice() { rand_stream.close(); }

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    result_type operator()()
    {
        result_type randomNumber = 0;
        rand_stream.read(reinterpret_cast<char*>(&randomNumber), sizeof(randomNumber));
        if (rand_stream.fail()) {
            throw std::system_error(errno, std::generic_category(), "Failed to read from /dev/random");
        }

        return randomNumber;
    }

private:
    std::ifstream rand_stream;
};
