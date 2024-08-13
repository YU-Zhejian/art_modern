#include "NotImplementedException.hh"

namespace labw {
namespace art_modern {

    const char* NotImplementedException::what() const noexcept
    {
        return "This method have not been implemented. Please contact software maintainer.";
    }
} // art_modern
} // labw