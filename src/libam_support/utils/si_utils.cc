#include "libam_support/utils/si_utils.hh"

#include "libam_support/Constants.hh"

#include <fmt/format.h>

#include <cstddef>
#include <cstdint>
#include <string>

namespace labw::art_modern {
std::string format_with_commas(const std::size_t number)
{
    std::string num_str = std::to_string(number);
    int64_t insertPosition = static_cast<int64_t>(num_str.length()) - 3;

    while (insertPosition > 0) {
        num_str.insert(insertPosition, ",");
        insertPosition -= 3;
    }

    return num_str;
}
std::string to_si(const double number, const int precision)
{
    std::size_t unitIndex = 0;
    auto sizeInUnit = static_cast<double>(number);

    while (sizeInUnit >= 1024 && unitIndex < SI_UNITS_LENGTH) {
        sizeInUnit /= 1024;
        ++unitIndex;
    }
    return fmt::format("{:.{}f}{}", sizeInUnit, precision, SI_UNITS[unitIndex]);
}
} // namespace labw::art_modern
