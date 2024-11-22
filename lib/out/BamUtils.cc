#include "BamUtils.hh"
#include "CExceptionsProxy.hh"
#include "art_modern_config.h"
#include "art_modern_constants.hh"
#include "utils/mpi_utils.hh"
#include "utils/seq_utils.hh"
#include <boost/log/trivial.hpp>

namespace labw::art_modern {
std::string BamUtils::generate_oa_tag(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar)
{
    const hts_pos_t pos = pwa.pos_on_contig + 1; // SAM is 1-based
    const auto strand = pwa.is_plus_strand ? '+' : '-';
    const auto cigar_str = cigar_arr_to_str(cigar);
    std::ostringstream oss;
    const auto nm_tag = "";
    oss << pwa.contig_name << ',' << pos << ',' << strand << ',' << cigar_str << ',' << MAPQ_MAX << ',' << nm_tag
        << ';';
    return oss.str();
}
BamUtils::BamUtils(const SamOptions& sam_options)
    : sam_options_(sam_options)
{
}
std::string BamUtils::generate_oa_tag(
    const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const int32_t nm_tag)
{
    const hts_pos_t pos = pwa.pos_on_contig + 1; // SAM is 1-based
    const auto strand = pwa.is_plus_strand ? '+' : '-';
    const auto cigar_str = cigar_arr_to_str(cigar);
    std::ostringstream oss;
    oss << pwa.contig_name << ',' << pos << ',' << strand << ',' << cigar_str << ',' << MAPQ_MAX << ','
        << std::to_string(nm_tag) << ';';
    return oss.str();
}
std::pair<int32_t, std::string> BamUtils::generate_nm_md_tag(
    const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar)
{
    hts_pos_t pos_on_query = 0;
    uint32_t matched = 0;
    hts_pos_t pos_on_ref = 0;
    std::stringstream md_str_ss;
    int32_t nm = 0;
    uint32_t this_cigar_len;
    uint32_t this_cigar_ops;

    for (auto i = 0; i < cigar.size(); ++i) {
        this_cigar_len = bam_cigar_oplen(cigar[i]);
        this_cigar_ops = bam_cigar_op(cigar[i]);
        if (this_cigar_ops == BAM_CEQUAL) {
            matched += this_cigar_len;
            pos_on_query += this_cigar_len;
            pos_on_ref += this_cigar_len;
        } else if (this_cigar_ops == BAM_CMATCH || this_cigar_ops == BAM_CDIFF) {
            for (int j = 0; j < this_cigar_len; ++j) {
                if (pwa.query[pos_on_query] == pwa.ref[pos_on_ref]) { // a match. TODO: Support IUB code
                    ++matched;
                } else {
                    md_str_ss << matched << static_cast<char>(std::toupper(pwa.ref[pos_on_ref]));
                    matched = 0;
                    ++nm;
                }
                pos_on_query++;
                pos_on_ref++;
            }
        } else if (this_cigar_ops == BAM_CDEL) {
            md_str_ss << matched << '^';
            for (int j = 0; j < this_cigar_len; ++j) {
                md_str_ss << static_cast<char>(std::toupper(pwa.ref[pos_on_ref]));
                pos_on_ref++;
                nm++;
            }
            matched = 0;
        } else if (this_cigar_ops == BAM_CINS || this_cigar_ops == BAM_CSOFT_CLIP) {
            pos_on_query += this_cigar_len;
            if (this_cigar_ops == BAM_CINS) {
                nm += static_cast<int32_t>(this_cigar_len);
            }
        } else if (this_cigar_ops == BAM_CREF_SKIP) {
            pos_on_ref += this_cigar_len;
        }
    }
    md_str_ss << matched;
    const auto md_str = md_str_ss.str();
    return { nm, md_str };
}
void assert_correct_cigar(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar)
{
#ifdef CEU_CM_IS_DEBUG
    const auto n_cigar = static_cast<int>(cigar.size());
    const auto cigar_qlen = bam_cigar2qlen(n_cigar, cigar.data());
    const auto cigar_rlen = bam_cigar2rlen(n_cigar, cigar.data());

    hts_pos_t pos_on_read = 0;
    hts_pos_t pos_on_ref = 0;
    uint32_t this_cigar_ops;
    uint32_t this_cigar_len;
    uint8_t this_cigar_type;

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
    bool is_equal;
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
    if (pos_on_read != pwa.query.length() || pos_on_ref != pwa.ref.length()) {
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
void BamTags::patch(bam1_t* record) const
{
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        CExceptionsProxy::assert_numeric(bam_aux_append(record, tag_name.c_str(), tag_type, tag_len, tag_data.get()),
            USED_HTSLIB_NAME, "Failed to add tag to read", false, CExceptionsProxy::EXPECTATION::ZERO);
    }
}
size_t BamTags::size() const
{
    size_t size = 0;
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        size += (tag_len + 3);
    }
    return size;
}
void BamTags::add_string(const std::string& key, const std::string& value)
{
    size_t len = value.size();
    data_type data(new uint8_t[len + 1]);
    std::copy(value.begin(), value.end(), data.get());
    data[static_cast<std::ptrdiff_t>(len)] = '\0';
    tags_.emplace_back(key, 'Z', len + 1, data);
}
void BamTags::add_int_c(const std::string& key, int8_t value)
{
    data_type data(new uint8_t[1]);
    data[0] = static_cast<uint8_t>(value);
    tags_.emplace_back(key, 'c', 1, data);
}
void BamTags::add_int_C(const std::string& key, uint8_t value)
{
    data_type data(new uint8_t[1]);
    data[0] = value;
    tags_.emplace_back(key, 'C', 1, data);
}
void BamTags::add_int_s(const std::string& key, int16_t value)
{
    data_type data(new uint8_t[2]);
    std::memcpy(data.get(), &value, 2);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 2);
#endif
    tags_.emplace_back(key, 's', 2, data);
}
void BamTags::add_int_S(const std::string& key, uint16_t value)
{
    data_type data(new uint8_t[2]);
    std::memcpy(data.get(), &value, 2);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 2);
#endif
    tags_.emplace_back(key, 'S', 2, data);
}
void BamTags::add_int_i(const std::string& key, int32_t value)
{
    data_type data(new uint8_t[4]);
    std::memcpy(data.get(), &value, 4);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 4);
#endif
    tags_.emplace_back(key, 'i', 4, data);
}
void BamTags::add_int_I(const std::string& key, uint32_t value)
{
    data_type data(new uint8_t[4]);
    std::memcpy(data.get(), &value, 4);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 4);
#endif
    tags_.emplace_back(key, 'I', 4, data);
}

}
