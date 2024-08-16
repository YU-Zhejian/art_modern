#pragma once

#include <iostream>

#include "out/BaseReadOutput.hh"
#include "stream/ThreadSafeFileStream.hh"
namespace labw {
namespace art_modern {
    class FastqReadOutput : public BaseReadOutput {
    public:
        explicit FastqReadOutput(const std::string& filename);
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

        void close() override;

        ~FastqReadOutput() override;

    private:
        ThreadSafeFileStream stream_;
    };

    class FastqReadOutputFactory : public BaseReadOutputFactory {
    public:
        void patch_options(boost::program_options::options_description& desc) override;
        std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
            const std::shared_ptr<BaseFastaFetch>& fasta_fetch) const override;
    };
}
}