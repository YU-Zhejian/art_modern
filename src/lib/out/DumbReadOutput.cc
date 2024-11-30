#include "DumbReadOutput.hh"

namespace labw::art_modern {
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
// labw