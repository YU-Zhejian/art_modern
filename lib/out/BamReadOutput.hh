#pragma once

#include "art_modern_constants.hh"
#include "fasta/BaseFastaFetch.hh"
#include "out/BaseReadOutput.hh"
#include <htslib/sam.h>
#include <memory>
#include <mutex>

namespace labw {
namespace art_modern {
    class SamReadOutputOptions {
    public:
        /**
         * Format version. Accepted format: `/^[0-9]+\.[0-9]+$`.
         */
        std::string HD_VN = "1.4";
        /**
         * Sorting order of alignments. Valid values: `unknown` (default), `unsorted`, `queryname` and `coordinate`.
         */
        std::string HD_SO = "unsorted";

        /**
         * Program record identifier.
         */
        std::string PG_ID = "01";
        /**
         * Program name.
         */
        std::string PG_PN = "art_modern";
        /**
         * Command line.
         */
        std::string PG_CL;
        /**
         * Program version.
         */
        std::string PG_VN = "ART-" ART_VERSION "-ART_MODERN-" ART_MODERN_VERSION;

        bool use_m = false;

        /**
         * If `false`, will write SAM instead.
         */
        bool write_bam = true;
    };

    class BamReadOutput : public BaseReadOutput {
    public:
        BamReadOutput(const std::string& filename, const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
            const SamReadOutputOptions& sam_options);
        void writeSE(const PairwiseAlignment& pwa) override;
        void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
        void close() override;
        ~BamReadOutput() override;

    private:
        samFile* sam_file_;
        sam_hdr_t* sam_header_;
        const SamReadOutputOptions& sam_options_;
        std::mutex mutex_;
        bool is_closed_ = false;
    };

    class BamReadOutputFactory : public BaseReadOutputFactory {
    public:
        void patch_options(boost::program_options::options_description& desc) override;
        std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
            std::shared_ptr<BaseFastaFetch>& fasta_fetch) const override;

    private:
        SamReadOutputOptions sam_options_;
    };
}
}
