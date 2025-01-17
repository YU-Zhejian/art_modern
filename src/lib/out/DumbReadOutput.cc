#include "out/DumbReadOutput.hh"

#include "ds/PairwiseAlignment.hh"

namespace labw::art_modern {
void DumbReadOutput::close()
{
    // Do nothing!
}
DumbReadOutput::~DumbReadOutput() { DumbReadOutput::close(); }
void DumbReadOutput::writePE([[maybe_unused]] const PairwiseAlignment&, [[maybe_unused]] const PairwiseAlignment&)
{
    // Do nothing!
}
void DumbReadOutput::writeSE([[maybe_unused]] const PairwiseAlignment&)
{
    // Do nothing!
}
DumbReadOutput::DumbReadOutput() = default;

} // namespace labw::art_modern
// labw