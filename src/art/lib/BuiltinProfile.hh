#pragma once

#include <absl/base/attributes.h>

#include <cstddef>
#include <string>

namespace labw::art_modern {

struct BuiltinProfile {
    std::string r1_profile;
    std::string r2_profile;
    BuiltinProfile(const unsigned char* profile1_encoded, std::size_t profile1_encoded_size,
        const unsigned char* profile2_encoded, std::size_t profile2_encoded_size);
} ABSL_ATTRIBUTE_FUNC_ALIGN((64));
} // namespace labw::art_modern
