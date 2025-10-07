

// TODO:
// 1. Add argument parsing
// 2. Support PE BAM/SAM

#include "libam_support/CExceptionsProxy.hh"


#include <boost/log/trivial.hpp>


#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdlib>
#include <iostream>

using namespace labw::art_modern;


static void view_sam(samFile *in)
{

    auto* hdr = CExceptionsProxy::assert_not_null( sam_hdr_read(in), "HTSLib", "Failed to read SAM header.");

    auto*  b =CExceptionsProxy::assert_not_null( bam_init1(), "HTSLib", "Failed to init BAM record.");
        int ret = 0;
        while ((ret = sam_read1(in, hdr, b)) >= 0) {
            std::cout << bam_get_qname(b) << ": ";
            const auto* qual = bam_get_qual(b);
            if (qual[0] == 0xff){
                std::cout << "*";
            } else {
                for (int i = 0; i < b->core.l_qseq; ++i) {
                    std::cout << (char)(qual[i] + 33);
                }
            }
            std::cout << std::endl;
        }

    sam_hdr_destroy(hdr);
    bam_destroy1(b);
}

int main() {
    htsFile *file = CExceptionsProxy::assert_not_null(hts_open("/home/yuzj/Documents/pbsim3_modern/explore/benchmark_other_simulators/data/soybean_SRR16074289.bam", "r"), "HTSLib", "Failed to open HTS file.");
    const auto* format = hts_get_format(file);
    if(format->category != sequence_data){
        BOOST_LOG_TRIVIAL(error) << "File is not a sequence data file.";
        hts_close(file);
        return EXIT_FAILURE;
    }
    view_sam(static_cast<samFile*>(file));

    // Iterate through all sequences and print qualities

    hts_close(file);
    return EXIT_SUCCESS;
}
