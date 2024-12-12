#include "OutputDispatcher.hh"
#include "BamReadOutput.hh"
#include "FastqReadOutput.hh"
#include "HeadlessBamReadOutput.hh"
#include "PwaReadOutput.hh"
#include <boost/log/trivial.hpp>

namespace labw::art_modern {

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
BaseReadOutput* OutputDispatcherFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
{
    const auto output_dispatcher = new OutputDispatcher();
    for (auto const& factory : factories_) {
        output_dispatcher->add(factory->create(vm, fasta_fetch, args));
    }
    BOOST_LOG_TRIVIAL(info) << "All writers added";
    return output_dispatcher;
}
void OutputDispatcherFactory::add(std::shared_ptr<BaseReadOutputFactory> factory)
{
    factories_.emplace_back(std::move(factory));
}

OutputDispatcherFactory::~OutputDispatcherFactory() = default;
OutputDispatcherFactory get_output_dispatcher_factory() noexcept
{
    OutputDispatcherFactory out_dispatcher_factory;
    out_dispatcher_factory.add(std::make_shared<PwaReadOutputFactory>());
    out_dispatcher_factory.add(std::make_shared<FastqReadOutputFactory>());
    out_dispatcher_factory.add(std::make_shared<BamReadOutputFactory>());
    out_dispatcher_factory.add(std::make_shared<HeadlessBamReadOutputFactory>());
    return out_dispatcher_factory;
}
}