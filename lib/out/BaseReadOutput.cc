#include "BaseReadOutput.hh"
#include "NotImplementedException.hh"

namespace labw {
namespace art_modern {
    void BaseReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        throw NotImplementedException();
    }
    void BaseReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        throw NotImplementedException();
    }
    BaseReadOutput::~BaseReadOutput() = default;
} // art_modern
} // labw