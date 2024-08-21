#include "BaseReadOutput.hh"

namespace labw {
namespace art_modern {


    BaseReadOutput::~BaseReadOutput() = default;

    BaseReadOutputFactory::~BaseReadOutputFactory() = default;

    void DumbReadOutput::close()
    {
        // Do nothing!
    }
    DumbReadOutput::~DumbReadOutput() { DumbReadOutput::close(); }
    void DumbReadOutput::writePE(const PairwiseAlignment&, const PairwiseAlignment&)
    {
        // Do nothing!
    }
    void DumbReadOutput::writeSE(const PairwiseAlignment&)
    {
        // Do nothing!
    }
    DumbReadOutput::DumbReadOutput() = default;

} // art_modern
} // labw