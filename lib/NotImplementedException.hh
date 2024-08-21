#pragma once
#include <exception>

namespace labw {
namespace art_modern {

    struct NotImplementedException : public std::exception {
        const char* what() const noexcept override;
    };

} // art_modern
} // labw
