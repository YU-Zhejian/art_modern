#pragma once

#include <string>

namespace labw::art_modern {

struct BuiltinProfile {
    std::string r1_profile;
    std::string r2_profile;
    BuiltinProfile(const char* profile1_encoded, const char* profile2_encoded);
};
} // namespace labw::art_modern
