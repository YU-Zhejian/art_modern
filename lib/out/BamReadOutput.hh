#pragma once

#include "BamUtils.hh"
#include "art_modern_constants.hh"
#include "fasta/BaseFastaFetch.hh"
#include "out/BaseReadOutput.hh"
#include <htslib/sam.h>

#include <mutex>

namespace labw {
namespace art_modern {

    class BamReadOutput : public BaseReadOutput {
    public:
        BamReadOutput(BamReadOutput&& other) = delete;
        BamReadOutput(const BamReadOutput&) = delete;
        BamReadOutput& operator=(BamReadOutput&&) = delete;
        BamReadOutput& operator=(const BamReadOutput&) = delete;

        BamReadOutput(const std::string& filename, const BaseFastaFetch* fasta_fetch, SamOptions sam_options);
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
        void close() override;
        ~BamReadOutput() override;

    private:
        samFile* sam_file_;
        sam_hdr_t* sam_header_;
        SamOptions sam_options_;
        std::mutex mutex_;
        bool is_closed_ = false;
    };

    class BamReadOutputFactory : public BaseReadOutputFactory {
    public:
        void patch_options(boost::program_options::options_description& desc) const override;
        BaseReadOutput* create(
            const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch) const override;
        ~BamReadOutputFactory() override;

    private:
        SamOptions sam_options_;
    };
}
}
