#include "BamUtils.hh"
#include "CExceptionsProxy.hh"
#include "art_modern_config.h"
#include "art_modern_constants.hh"
#include "utils/mpi_utils.hh"
#include "utils/seq_utils.hh"
#include <boost/log/trivial.hpp>

namespace labw::art_modern {
std::string BamUtils::generate_oa_tag(const PairwiseAlignment& pwa) const
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
BamUtils::BamUtils(const SamOptions& sam_options)
    : sam_options_(sam_options)
{
}

void assert_correct_cigar(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const uint32_t* cigar_c_arr)
{
#ifdef CEU_CM_IS_DEBUG
    auto cigar_qlen = bam_cigar2qlen(static_cast<int>(cigar.size()), cigar_c_arr);
    auto cigar_rlen = bam_cigar2rlen(static_cast<int>(cigar.size()), cigar_c_arr);

    if (cigar_qlen != pwa.query.length()) {
        BOOST_LOG_TRIVIAL(error) << "Cigar length mismatch with query: " << cigar_qlen << " != " << pwa.query.length();
        goto err;
    }

    if (cigar_rlen != pwa.ref.length()) {
        BOOST_LOG_TRIVIAL(error) << "Cigar length mismatch with ref: " << cigar_rlen << " != " << pwa.ref.length();
        goto err;
    }
    if (pwa.query.length() != pwa.qual.length()) {
        BOOST_LOG_TRIVIAL(error) << "Qual length mismatch with query: " << pwa.qual.length()
                                 << " != " << pwa.query.length();
        goto err;
    }
    return;
err:
    BOOST_LOG_TRIVIAL(error) << "Cigar  : " << cigar_arr_to_str(cigar) << " (QLEN=" << cigar_qlen
                             << ", RLEN=" << cigar_rlen << ")";
    BOOST_LOG_TRIVIAL(error) << "Query  : " << pwa.query;
    BOOST_LOG_TRIVIAL(error) << "Qual   : " << pwa.qual;
    BOOST_LOG_TRIVIAL(error) << "Ref    : " << pwa.ref;
    BOOST_LOG_TRIVIAL(error) << "AQuery : " << pwa.aligned_query;
    BOOST_LOG_TRIVIAL(error) << "ARef   : " << pwa.aligned_ref;
    abort_mpi();
#endif
}
void fill_md_nm_tag(bam1_t* b, const PairwiseAlignment& pwa)
{
    const auto cigar_encoded = bam_get_cigar(b);
    const auto c = &b->core;
    hts_pos_t qpos = 0;
    int matched = 0;
    hts_pos_t rpos = 0;
    std::stringstream md_str_ss;
    int32_t nm = 0;

    for (int i = 0; i < c->n_cigar; ++i) {
        int oplen = bam_cigar_oplen(cigar_encoded[i]);
        int op = bam_cigar_op(cigar_encoded[i]);
        char c1;
        char c2;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < oplen; ++j) {
                c1 = pwa.query[qpos];
                c2 = pwa.ref[rpos];
                if (c1 == c2) { // a match. TODO: Support IUB code
                    ++matched;
                } else {
                    md_str_ss << matched << static_cast<char>(std::toupper(pwa.ref[rpos]));
                    matched = 0;
                    ++nm;
                }
                qpos++;
                rpos++;
            }
        } else if (op == BAM_CDEL) {
            md_str_ss << matched << '^';
            for (int j = 0; j < oplen; ++j) {
                md_str_ss << static_cast<char>(std::toupper(pwa.ref[rpos]));
                rpos++;
                nm++;
            }
            matched = 0;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qpos += oplen;
            if (op == BAM_CINS) {
                nm += oplen;
            }
        } else if (op == BAM_CREF_SKIP) {
            rpos += oplen;
        }
    }
    md_str_ss << matched;

    // update NM
    CExceptionsProxy::assert_numeric(bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm), USED_HTSLIB_NAME,
        "Failed to add NM tag", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    // update MD
    const auto md_str = md_str_ss.str();
    const auto md_cstr = md_str.c_str();
    CExceptionsProxy::assert_numeric(
        bam_aux_append(b, "MD", 'Z', static_cast<int>(md_str.size() + 1), (uint8_t*)md_cstr), USED_HTSLIB_NAME,
        "Failed to add NM tag", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
}
}
