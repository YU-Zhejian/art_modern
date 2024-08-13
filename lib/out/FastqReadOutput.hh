#pragma once

#include <iostream>

#include "out/BaseReadOutput.hh"
#include "stream/ThreadSafeFileStream.hh"
namespace labw {
namespace art_modern {

    /**
     * This class will **NOT** close underlying stream!
     */
    class FastqReadOutput : public BaseReadOutput {
    public:
        explicit FastqReadOutput(ThreadSafeFileStream& stream);
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
        ~FastqReadOutput() override;

    private:
        ThreadSafeFileStream& stream_;
    };

}
}