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

    void OutputDispatcher::add(const std::shared_ptr<BaseReadOutput>& output) { outputs_.emplace_back(output); }

    void OutputDispatcherFactory::patch_options(boost::program_options::options_description& desc)
    {
        for (auto const& factory : factories_) {
            factory->patch_options(desc);
        }
    }
    std::shared_ptr<BaseReadOutput> OutputDispatcherFactory::create(
        const boost::program_options::variables_map& vm, const std::shared_ptr<BaseFastaFetch>& fasta_fetch) const
    {
        auto output_dispatcher = std::make_shared<OutputDispatcher>();
        for (auto const& factory : factories_) {
            output_dispatcher->add(factory->create(vm, fasta_fetch));
        }
        return output_dispatcher;
    }
    void OutputDispatcherFactory::add(const std::shared_ptr<BaseReadOutputFactory>& factory)
    {
        factories_.emplace_back(factory);
    }

}
}