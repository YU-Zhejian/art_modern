#pragma once
#include <memory>
#include <vector>

#include "BaseReadOutput.hh"

namespace labw {
namespace art_modern {

    class OutputDispatcher : public BaseReadOutput {

    public:
        void add(std::shared_ptr<BaseReadOutput> output);

        void writeSE(const PairwiseAlignment& pwa) override;

        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

        void close() override;

        ~OutputDispatcher() override;

    private:
        std::vector<std::shared_ptr<BaseReadOutput>> outputs_;
    };

}
}