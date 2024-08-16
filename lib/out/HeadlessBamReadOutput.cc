

#include "CExceptionsProxy.hh"



#include <boost/algorithm/string/join.hpp>

#include <boost/log/trivial.hpp>



#include "HeadlessBamReadOutput.hh"



#include <seq_utils.hh>



namespace labw {

namespace art_modern {

    HeadlessBamReadOutput::HeadlessBamReadOutput(const std::string& filename, const SamReadOutputOptions& sam_options)

        : sam_options_(sam_options)

    {

        std::unique_lock<std::mutex> rhs_lk(mutex_);



        sam_file_ = (samFile*)CExceptionsProxy::requires_not_null(

            sam_open(filename.c_str(), sam_options_.write_bam ? "wb" : "w"), USED_HTSLIB_NAME,

            "Failed to open SAM file");

        sam_header_ = (sam_hdr_t*)CExceptionsProxy::requires_not_null(

            sam_hdr_init(), USED_HTSLIB_NAME, "Faield to initialize SAM header");

        CExceptionsProxy::requires_numeric(sam_hdr_add_line(sam_header_, "HD", "VN", sam_options_.HD_VN.c_str(), "SO",

                                               sam_options_.HD_SO.c_str(), NULL),

            USED_HTSLIB_NAME, "Failed to add HD header line");

        CExceptionsProxy::requires_numeric(sam_hdr_add_line(sam_header_, "PG", "ID", sam_options_.PG_ID.c_str(), "PN",

                                               sam_options_.PG_PN.c_str(), "CL", sam_options_.PG_CL.c_str(), NULL),

            USED_HTSLIB_NAME, "Failed to add PG header line", false, CExceptionsProxy::EXPECTATION::ZERO);



        CExceptionsProxy::requires_numeric(

            sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");

        BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";

    }

    void HeadlessBamReadOutput::writeSE(const PairwiseAlignment& pwa)

    {

        if (is_closed_) {

            return;

        }

        std::unique_lock<std::mutex> rhs_lk(mutex_);



        auto sam_record = (bam1_t*)CExceptionsProxy::requires_not_null(

            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");

        auto rlen = static_cast<long>(pwa.query.size());



        auto seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);

        auto qual = pwa.qual;

        if (!pwa.is_plus_strand) {

            std::reverse(qual.begin(), qual.end());

        }

        auto oa_tag = generate_oa_tag(pwa);

        auto tag_len = 2 + 1 + oa_tag.size() + 1;

        CExceptionsProxy::requires_numeric(bam_set1(sam_record,

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               0, // Alignment info moved to OA tag

                                               static_cast<uint32_t*>(nullptr), // Alignment info moved to OA tag

                                               0, // Unset for SE reads

                                               0, // Unset for SE reads

                                               0, // Unset for SE reads

                                               rlen, seq.c_str(), qual.c_str(), tag_len),

            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);



        bam_aux_update_str(sam_record, "OA", oa_tag.length(), oa_tag.c_str());

        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record), USED_HTSLIB_NAME,

            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    }

    void HeadlessBamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)

    {

        // TODO

    }



    void HeadlessBamReadOutput::close()

    {

        if (is_closed_) {

            return;

        }

        std::unique_lock<std::mutex> rhs_lk(mutex_);

        sam_close(sam_file_);

        is_closed_ = true;

    }



    HeadlessBamReadOutput::~HeadlessBamReadOutput() { HeadlessBamReadOutput::close(); }

    std::string HeadlessBamReadOutput::generate_oa_tag(const PairwiseAlignment& pwa) const

    {

        auto cigar = pwa.generate_cigar_array(sam_options_.use_m);



        auto seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);

        hts_pos_t pos = pwa.align_contig_start + 1; // SAM is 1-based

        auto strand = pwa.is_plus_strand ? '+' : '-';



        auto cigar_str = cigar_arr_to_str(cigar);

        std::ostringstream oss;

        auto nm_tag = "";

        oss << pwa.contig_name << ',' << pos << ',' << strand << ',' << cigar_str << ',' << MAPQ_MAX << ',' << nm_tag

            << ';';

        return oss.str();

    }

}

}

