#pragma once
#include "PairwiseAlignment.hh"
#include "fasta/BaseFastaFetch.hh"
#include <boost/program_options.hpp>

namespace labw {
namespace art_modern {

    class BaseReadOutput {
    public:
        virtual void writeSE(const PairwiseAlignment& pwa);
        virtual void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2);
        virtual void close();
        virtual ~BaseReadOutput();
    };

    class BaseReadOutputFactory {
    public:
        virtual void patch_options(boost::program_options::options_description& desc);
        virtual std::shared_ptr<BaseReadOutput> create(
            const boost::program_options::variables_map& vm, std::shared_ptr<BaseFastaFetch>& fasta_fetch) const;
    };

    class DumbReadOutput final : public BaseReadOutput {
    public:
        DumbReadOutput();
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
        void close() override;
        ~DumbReadOutput() override;
    };

} // art_modern
} // labw
