#include "art_modern_config.h"

#include "libam/bam/BamUtils.hh"

#include "libam/CExceptionsProxy.hh"
#include "libam/Constants.hh"
#include "libam/Dtypes.hh"
#include "libam/bam/BamOptions.hh"
#include "libam/bam/BamTypes.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/utils/mpi_utils.hh" // NOLINT
#include "libam/utils/seq_utils.hh"

#include <fmt/core.h>

#include <boost/log/trivial.hpp> // Used in DEBUG build of assert_correct_cigar

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

std::string BamUtils::generate_oa_tag(
    const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar, const int32_t nm_tag)
{
    return fmt::format("{},{},{},{},{},{};", pwa.contig_name, pwa.pos_on_contig + 1, pwa.is_plus_strand ? '+' : '-', cigar_arr_to_str(cigar), MAPQ_MAX, nm_tag);
}
std::pair<int32_t, std::string> BamUtils::generate_nm_md_tag(
    const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar)
{
    hts_pos_t pos_on_query = 0;
    uint32_t matched = 0;
    hts_pos_t pos_on_ref = 0;
    std::ostringstream md_str_ss;
    int32_t nm = 0;
    am_cigar_t this_cigar_len = 0;
    am_cigar_t this_cigar_ops = 0;

    for (const auto& this_cigar : cigar) {
        this_cigar_len = bam_cigar_oplen(this_cigar);
        this_cigar_ops = bam_cigar_op(this_cigar);
        if (this_cigar_ops == BAM_CEQUAL) {
            matched += this_cigar_len;
            pos_on_query += this_cigar_len;
            pos_on_ref += this_cigar_len;
        } else if (this_cigar_ops == BAM_CDIFF) {
            for (decltype(this_cigar_len) j = 0; j < this_cigar_len; ++j) {
                md_str_ss << std::to_string(matched) << static_cast<char>(std::toupper(pwa.ref[pos_on_ref]));
                matched = 0;
                ++nm;
                pos_on_query++;
                pos_on_ref++;
            }
        } else if (this_cigar_ops == BAM_CMATCH) {
            for (decltype(this_cigar_len) j = 0; j < this_cigar_len; ++j) {
                if (pwa.query[pos_on_query] == pwa.ref[pos_on_ref]) {
                    ++matched;
                } else {
                    md_str_ss << std::to_string(matched) << static_cast<char>(std::toupper(pwa.ref[pos_on_ref]));
                    matched = 0;
                    ++nm;
                }
                pos_on_query++;
                pos_on_ref++;
            }
        } else if (this_cigar_ops == BAM_CDEL) {
            md_str_ss << std::to_string(matched) << '^' << pwa.ref.substr(pos_on_ref, this_cigar_len);
            pos_on_ref += this_cigar_len;
            nm += static_cast<int32_t>(this_cigar_len);
            matched = 0;
        } else if (this_cigar_ops == BAM_CINS) {
            pos_on_query += this_cigar_len;
            nm += static_cast<int32_t>(this_cigar_len);
        } else if (this_cigar_ops == BAM_CSOFT_CLIP) {
            pos_on_query += this_cigar_len;
        } else if (this_cigar_ops == BAM_CREF_SKIP) {
            pos_on_ref += this_cigar_len;
        }
    }
    md_str_ss << std::to_string(matched);
    const auto md_str = md_str_ss.str();

    return { nm, md_str };
}

sam_hdr_t* BamUtils::init_header(const BamOptions& sam_options)
{
    auto* const sam_header
        = CExceptionsProxy::assert_not_null(sam_hdr_init(), USED_HTSLIB_NAME, "Faield to initialize SAM header");

    CExceptionsProxy::assert_numeric(
        sam_hdr_add_line(sam_header, "HD", "VN", sam_options.HD_VN.c_str(), "SO", sam_options.HD_SO.c_str(), NULL),
        USED_HTSLIB_NAME, "Failed to add HD header line");
    CExceptionsProxy::assert_numeric(sam_hdr_add_line(sam_header, "PG", "ID", sam_options.PG_ID.c_str(), "PN",
                                         sam_options.PG_PN.c_str(), "CL", sam_options.PG_CL.c_str(), NULL),
        USED_HTSLIB_NAME, "Failed to add PG header line", false, CExceptionsProxy::EXPECTATION::ZERO);
    return sam_header;
}
samFile* BamUtils::open_file(const std::string& filename, const BamOptions& sam_options)
{
    std::string mode;
    if (sam_options.write_bam) {
        if (!ends_with(filename, ".bam")) {
            BOOST_LOG_TRIVIAL(warning) << "BAM file name was not end with .bam: " << filename;
        }
        mode += "wb";
        mode += std::to_string(sam_options.compress_level);
    } else {
        if (!ends_with(filename, ".sam")) {
            BOOST_LOG_TRIVIAL(warning) << "SAM file name was not end with .bam: " << filename;
        }
        mode += "wh";
    }

    auto* const retv = CExceptionsProxy::assert_not_null(
        sam_open(filename.c_str(), mode.c_str()), USED_HTSLIB_NAME, "Failed to open SAM file");
    CExceptionsProxy::assert_numeric(hts_set_threads(retv, sam_options.hts_io_threads), USED_HTSLIB_NAME,
        "Failed to set writer thread number", false, CExceptionsProxy::EXPECTATION::ZERO);
    return retv;
}

bam1_t_uptr BamUtils::init_uptr()
{
    return bam1_t_uptr { CExceptionsProxy::assert_not_null(
        bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record") };
}

void BamUtils::assert_correct_cigar(
    [[maybe_unused]] const PairwiseAlignment& pwa, [[maybe_unused]] const std::vector<am_cigar_t>& cigar)
{
#ifdef CEU_CM_IS_DEBUG
    const auto n_cigar = static_cast<int>(cigar.size());
    const auto cigar_qlen = bam_cigar2qlen(n_cigar, cigar.data());
    const auto cigar_rlen = bam_cigar2rlen(n_cigar, cigar.data());

    hts_pos_t pos_on_read = 0;
    hts_pos_t pos_on_ref = 0;
    am_cigar_t this_cigar_ops = 0;
    am_cigar_t this_cigar_len = 0;
    am_cigar_type_t this_cigar_type = 0;
    bool is_equal = false;

    const auto qlen = static_cast<hts_pos_t>(pwa.query.length());
    const auto rlen = static_cast<hts_pos_t>(pwa.ref.length());
    const auto qual_len = static_cast<hts_pos_t>(pwa.qual.length());

    if (cigar_qlen != qlen) {
        BOOST_LOG_TRIVIAL(error) << "Cigar length mismatch with query: " << cigar_qlen << " != " << qlen;
        goto err;
    }
    if (cigar_rlen != rlen) {
        BOOST_LOG_TRIVIAL(error) << "Cigar length mismatch with ref: " << cigar_rlen << " != " << rlen;
        goto err;
    }
    if (qlen != qual_len) {
        BOOST_LOG_TRIVIAL(error) << "Qual length mismatch with query: " << qual_len << " != " << qual_len;
        goto err;
    }
    for (auto cigar_idx = 0; cigar_idx < n_cigar; ++cigar_idx) {
        this_cigar_ops = bam_cigar_op(cigar[cigar_idx]);
        this_cigar_len = bam_cigar_oplen(cigar[cigar_idx]);
        this_cigar_type = bam_cigar_type(this_cigar_ops);

        switch (this_cigar_type) {
        case CONSUME_QUERY_AND_REFERENCE:
            if (this_cigar_ops == BAM_CEQUAL || this_cigar_ops == BAM_CDIFF) {
                is_equal
                    = std::strncmp(pwa.query.c_str() + pos_on_read, pwa.ref.c_str() + pos_on_ref, this_cigar_len) == 0;
                if ((this_cigar_ops == BAM_CEQUAL) != is_equal) {
                    BOOST_LOG_TRIVIAL(error) << "Query/ref X/= error in CIGAR. QPOS=" << pos_on_read
                                             << " RPOS=" << pos_on_ref << " CIGAR_ID=" << cigar_idx << " ("
                                             << this_cigar_len << bam_cigar_opchr(this_cigar_ops) << ")";
                    goto err;
                }
            }
            pos_on_read += this_cigar_len;
            pos_on_ref += this_cigar_len;
            break;
        case CONSUME_NEITHER_QUERY_NOR_REFERENCE:
            break;
        case CONSUME_QUERY:
            pos_on_read += this_cigar_len;
            break;
        case CONSUME_REFERENCE:
            pos_on_ref += this_cigar_len;
            break;
        default: // Error!
            abort_mpi();
        }
    }
    if (pos_on_read != qlen || pos_on_ref != rlen) {
        BOOST_LOG_TRIVIAL(error) << "Cigar length mismatch with query and ref: " << pos_on_read
                                 << " != " << pwa.query.length() << " or " << pos_on_ref << " != " << pwa.ref.length();
        goto err;
    }
    return;
err:
    BOOST_LOG_TRIVIAL(error) << "Read   : " << pwa.read_name;
    BOOST_LOG_TRIVIAL(error) << "Contig : " << pwa.contig_name << ":" << pwa.pos_on_contig << ":"
                             << (pwa.is_plus_strand ? "+" : "-");

    BOOST_LOG_TRIVIAL(error) << "Cigar  : " << cigar_arr_to_str(cigar) << " (QLEN=" << cigar_qlen
                             << ", RLEN=" << cigar_rlen << ")";
    BOOST_LOG_TRIVIAL(error) << "Query  : " << pwa.query;
    BOOST_LOG_TRIVIAL(error) << "Ref    : " << pwa.ref;
    BOOST_LOG_TRIVIAL(error) << "Qual   : " << pwa.qual;
    BOOST_LOG_TRIVIAL(error) << "AQuery : " << pwa.aligned_query;
    BOOST_LOG_TRIVIAL(error) << "ARef   : " << pwa.aligned_ref;
    abort_mpi();
#endif
}

} // namespace labw::art_modern
