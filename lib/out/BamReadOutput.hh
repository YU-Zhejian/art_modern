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
        BamReadOutput(
            const std::string& filename, BaseFastaFetch* fasta_fetch, const SamReadOutputOptions& sam_options);
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
        void close() override;
        ~BamReadOutput() override;

    private:
        samFile* sam_file_;
        sam_hdr_t* sam_header_;
        SamReadOutputOptions sam_options_;
        std::mutex mutex_;
        bool is_closed_ = false;
        BamUtils bam_utils_;
    };

    class BamReadOutputFactory : public BaseReadOutputFactory {
    public:
        void patch_options(boost::program_options::options_description& desc) override;
        BaseReadOutput* create(
            const boost::program_options::variables_map& vm, BaseFastaFetch* fasta_fetch) const override;
        ~BamReadOutputFactory();

    private:
        SamReadOutputOptions sam_options_;
    };
}
}
