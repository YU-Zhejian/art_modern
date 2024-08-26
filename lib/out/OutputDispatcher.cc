#include "OutputDispatcher.hh"
#include <boost/log/trivial.hpp>

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
    }

    OutputDispatcher::~OutputDispatcher()
    {
        OutputDispatcher::close();

        for (const auto& output : outputs_) {
            delete output;
        }
    }

    void OutputDispatcher::add(BaseReadOutput* output) { outputs_.emplace_back(output); }

    void OutputDispatcherFactory::patch_options(boost::program_options::options_description& desc) const
    {
        for (auto const& factory : factories_) {
            factory->patch_options(desc);
        }
    }
    BaseReadOutput* OutputDispatcherFactory::create(
        const boost::program_options::variables_map& vm, const BaseFastaFetch *fasta_fetch) const
    {
        auto output_dispatcher = new OutputDispatcher();
        for (auto const& factory : factories_) {
            output_dispatcher->add(factory->create(vm, fasta_fetch));
        }
        BOOST_LOG_TRIVIAL(info) << "All writers added";
        return output_dispatcher;
    }
    void OutputDispatcherFactory::add(BaseReadOutputFactory* factory) { factories_.emplace_back(factory); }

    OutputDispatcherFactory::~OutputDispatcherFactory()
    {
        for (auto& factory : factories_) {
            delete factory;
        }
    }

}
}