#pragma once

#include <absl/base/attributes.h>

#include <string>

namespace labw::art_modern {

struct BuiltinProfile {
    std::string r1_profile;
    std::string r2_profile;
    BuiltinProfile(const char* profile1_encoded, const char* profile2_encoded);
} ABSL_ATTRIBUTE_FUNC_ALIGN((64));
} // namespace labw::art_modern
