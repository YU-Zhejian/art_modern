#pragma once
#include "PairwiseAlignment.hh"

namespace labw {
namespace art_modern {

    class BaseReadOutput {
    public:
        virtual void writeSE(const PairwiseAlignment& pwa);
        virtual void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2);
        virtual void close();
        virtual ~BaseReadOutput();
    };

} // art_modern
} // labw
