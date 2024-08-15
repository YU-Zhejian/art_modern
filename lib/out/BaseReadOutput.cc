#include "BaseReadOutput.hh"
#include "NotImplementedException.hh"

namespace labw {
namespace art_modern {
    void BaseReadOutput::writeSE(const PairwiseAlignment& pwa) { throw NotImplementedException(); }
    void BaseReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        throw NotImplementedException();
    }

    void BaseReadOutput::close() { throw NotImplementedException(); }

    BaseReadOutput::~BaseReadOutput() = default;
    void BaseReadOutputFactory::patch_options(boost::program_options::options_description& desc)
    {
        throw NotImplementedException();
    }
    std::shared_ptr<BaseReadOutput> BaseReadOutputFactory::create(
        const boost::program_options::variables_map& vm, std::shared_ptr<BaseFastaFetch>& fasta_fetch) const
    {
        throw NotImplementedException();
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