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
std::string BamUtils::generate_oa_tag(
    const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const int32_t nm_tag)
{
    hts_pos_t pos = pwa.align_contig_start + 1; // SAM is 1-based
    auto strand = pwa.is_plus_strand ? '+' : '-';
    auto cigar_str = cigar_arr_to_str(cigar);
    std::ostringstream oss;
    oss << pwa.contig_name << ',' << pos << ',' << strand << ',' << cigar_str << ',' << MAPQ_MAX << ','
        << std::to_string(nm_tag) << ';';
    return oss.str();
}
std::pair<int32_t, std::string> BamUtils::generate_nm_md_tag(
    const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar)
{
    hts_pos_t qpos = 0;
    int matched = 0;
    hts_pos_t rpos = 0;
    std::stringstream md_str_ss;
    int32_t nm = 0;
    uint32_t oplen;
    uint32_t op;

    for (auto i = 0; i < cigar.size(); ++i) {
        oplen = bam_cigar_oplen(cigar[i]);
        op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < oplen; ++j) {
                if (pwa.query[qpos] == pwa.ref[rpos]) { // a match. TODO: Support IUB code
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
    const auto md_str = md_str_ss.str();
    return { nm, md_str };
}
void assert_correct_cigar(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar)
{
#ifdef CEU_CM_IS_DEBUG
    auto cigar_qlen = bam_cigar2qlen(static_cast<int>(cigar.size()), cigar.data());
    auto cigar_rlen = bam_cigar2rlen(static_cast<int>(cigar.size()), cigar.data());

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
    std::shared_ptr<uint8_t[]> data(new uint8_t[len + 1]);
    std::copy(value.begin(), value.end(), data.get());
    data[static_cast<std::ptrdiff_t>(len)] = '\0';
    tags_.emplace_back(key, 'Z', len + 1, data);
}
void BamTags::add_int_c(const std::string& key, int8_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[1]);
    data[0] = static_cast<uint8_t>(value);
    tags_.emplace_back(key, 'c', 1, data);
}
void BamTags::add_int_C(const std::string& key, uint8_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[1]);
    data[0] = value;
    tags_.emplace_back(key, 'C', 1, data);
}
void BamTags::add_int_s(const std::string& key, int16_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[2]);
    std::memcpy(data.get(), &value, 2);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 2);
#endif
    tags_.emplace_back(key, 's', 2, data);
}
void BamTags::add_int_S(const std::string& key, uint16_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[2]);
    std::memcpy(data.get(), &value, 2);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 2);
#endif
    tags_.emplace_back(key, 'S', 2, data);
}
void BamTags::add_int_i(const std::string& key, int32_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[4]);
    std::memcpy(data.get(), &value, 4);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 4);
#endif
    tags_.emplace_back(key, 'i', 4, data);
}
void BamTags::add_int_I(const std::string& key, uint32_t value)
{
    std::shared_ptr<uint8_t[]> data(new uint8_t[4]);
    std::memcpy(data.get(), &value, 4);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 4);
#endif
    tags_.emplace_back(key, 'I', 4, data);
}

}
