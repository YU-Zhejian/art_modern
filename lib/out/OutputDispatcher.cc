#include "OutputDispatcher.hh"
namespace labw {
namespace art_modern {

    void OutputDispatcher::writeSE(const PairwiseAlignment& pwa)
    {
        for (const auto& output : outputs_) {
            output->writeSE(pwa);
        }
    }

    void OutputDispatcher::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        for (const auto& output : outputs_) {
            output->writePE(pwa1, pwa2);
        }
    }

    void OutputDispatcher::close()
    {
        for (const auto& output : outputs_) {
            output->close();
        }
        outputs_.clear();
    }

    OutputDispatcher::~OutputDispatcher() { OutputDispatcher::close(); }

    void OutputDispatcher::add(std::shared_ptr<BaseReadOutput> output) { outputs_.emplace_back(output); }
}
}