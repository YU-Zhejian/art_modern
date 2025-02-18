#include "art/lib/BuiltinProfile.hh"

// #include <base64.h>

#include <cstring>

namespace labw::art_modern {
BuiltinProfile::BuiltinProfile(const char* profile1_encoded, const char* profile2_encoded)
    : r1_profile(std::string(profile1_encoded, std::strlen(profile1_encoded)))
    , r2_profile(profile2_encoded[0] != '\0' ? std::string(profile2_encoded, std::strlen(profile2_encoded)) : "")
{
}
} // namespace labw::art_modern
