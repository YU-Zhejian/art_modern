#pragma once
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/utils/class_macros_utils.hh"

namespace labw::art_modern {

class DumbReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(DumbReadOutput)
    DELETE_MOVE(DumbReadOutput)

    DumbReadOutput();
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;

    bool require_alignment() const override;

    ~DumbReadOutput() override;
};

} // namespace labw::art_modern
