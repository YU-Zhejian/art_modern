#pragma once

#include <htslib/hts.h>

#include <stdexcept>
namespace labw::art_modern {
class TidNotFound : public std::runtime_error {
public:
    TidNotFound()
        : std::runtime_error("Tid not found")
    {
    }
};
} // namespace labw::art_modern
