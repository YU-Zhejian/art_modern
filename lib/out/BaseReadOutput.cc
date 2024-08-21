#include "BaseReadOutput.hh"
#include "ExceptionUtils.hh"
#include "NotImplementedException.hh"
#include <boost/log/trivial.hpp>

namespace labw {
namespace art_modern {
    void BaseReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseReadOutput::writeSE?";
        throw_with_trace(NotImplementedException());
    }
    void BaseReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseReadOutput::writePE?";
        throw_with_trace(NotImplementedException());
    }

    void BaseReadOutput::close()
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseReadOutput::close?";
        throw_with_trace(NotImplementedException());
    }

    BaseReadOutput::~BaseReadOutput() = default;
    void BaseReadOutputFactory::patch_options(boost::program_options::options_description& desc)
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseReadOutputFactory::patch_options?";
        throw_with_trace(NotImplementedException());
    }
    std::shared_ptr<BaseReadOutput> BaseReadOutputFactory::create(
        const boost::program_options::variables_map& vm, const std::shared_ptr<BaseFastaFetch>& fasta_fetch) const
    {
        BOOST_LOG_TRIVIAL(error) << "Some method calls BaseReadOutputFactory::create?";
        throw_with_trace(NotImplementedException());
    }
    void DumbReadOutput::close()
    {
        // Do nothing!
    }
    DumbReadOutput::~DumbReadOutput() { DumbReadOutput::close(); }
    void DumbReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        // Do nothing!
    }
    void DumbReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        // Do nothing!
    }
    DumbReadOutput::DumbReadOutput() = default;

} // art_modern
} // labw