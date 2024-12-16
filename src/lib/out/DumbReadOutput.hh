#pragma once
#include "PairwiseAlignment.hh"
#include "out/BaseReadOutput.hh"

namespace labw::art_modern {

class DumbReadOutput final : public BaseReadOutput {
public:
    DumbReadOutput(DumbReadOutput&& other) = delete;
    DumbReadOutput(const DumbReadOutput&) = delete;
    DumbReadOutput& operator=(DumbReadOutput&&) = delete;
    DumbReadOutput& operator=(const DumbReadOutput&) = delete;

    DumbReadOutput();
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~DumbReadOutput() override;
};

} // art_modern
// labw
