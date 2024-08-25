#pragma once
#include "PairwiseAlignment.hh"
#include "fasta/BaseFastaFetch.hh"
#include <boost/program_options.hpp>

namespace labw {
namespace art_modern {

    class BaseReadOutput {
    public:
        BaseReadOutput(BaseReadOutput&& other) = delete;
        BaseReadOutput(const BaseReadOutput&) = delete;
        BaseReadOutput& operator=(BaseReadOutput&&) = delete;
        BaseReadOutput() = default;
        virtual void writeSE(const PairwiseAlignment& pwa) = 0;
        virtual void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) = 0;
        virtual void close() = 0;
        virtual ~BaseReadOutput();
    };

    class BaseReadOutputFactory {
    public:
        virtual void patch_options(boost::program_options::options_description& desc) const = 0;
        virtual BaseReadOutput* create(
            const boost::program_options::variables_map& vm, BaseFastaFetch* fasta_fetch) const
            = 0;
        virtual ~BaseReadOutputFactory();
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
