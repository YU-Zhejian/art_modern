#pragma once
#include "ds/PairwiseAlignment.hh"
#include "out/BaseReadOutput.hh"
#include "utils/class_macros_utils.hh"

namespace labw::art_modern {

class DumbReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(DumbReadOutput)
    DELETE_MOVE(DumbReadOutput)

    DumbReadOutput();
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~DumbReadOutput() override;
};

} // namespace labw::art_modern
