#pragma once
#include <memory>
#include <vector>

#include "BaseReadOutput.hh"

namespace labw {
namespace art_modern {

    class OutputDispatcher : public BaseReadOutput {

    public:
        void add(const std::shared_ptr<BaseReadOutput>& output);

        void writeSE(const PairwiseAlignment& pwa) override;

        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

        void close() override;

        ~OutputDispatcher() override;

    private:
        std::vector<std::shared_ptr<BaseReadOutput>> outputs_;
    };

    class OutputDispatcherFactory : public BaseReadOutputFactory {
    public:
        void add(const std::shared_ptr<BaseReadOutputFactory>& factory);
        void patch_options(boost::program_options::options_description& desc) override;
        std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
            std::shared_ptr<BaseFastaFetch>& fasta_fetch) const override;

    private:
        std::vector<std::shared_ptr<BaseReadOutputFactory>> factories_;
    };

}
}