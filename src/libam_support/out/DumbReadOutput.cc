#include "libam_support/out/DumbReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"

namespace labw::art_modern {
void DumbReadOutput::close()
{
    // Do nothing!
}
DumbReadOutput::~DumbReadOutput() { DumbReadOutput::close(); }
void DumbReadOutput::writePE(
    [[maybe_unused]] const PairwiseAlignment& /*pwa1*/, [[maybe_unused]] const PairwiseAlignment& /*pwa2*/)
{
    // Do nothing!
}
void DumbReadOutput::writeSE([[maybe_unused]] const PairwiseAlignment& /*pwa*/)
{
    // Do nothing!
}

bool DumbReadOutput::require_alignment() const { return false; }

DumbReadOutput::DumbReadOutput() = default;

} // namespace labw::art_modern
