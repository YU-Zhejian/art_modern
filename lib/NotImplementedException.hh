#pragma once
#include <exception>

namespace labw {
namespace art_modern {

    struct NotImplementedException : public std::exception {
        const char* what() const _GLIBCXX_TXN_SAFE_DYN _GLIBCXX_NOTHROW override;
    };

} // art_modern
} // labw
