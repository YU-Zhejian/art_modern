/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "art/lib/ArtRead.hh"

#include "art/lib/Rprob.hh"

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG
#include "libam_support/Constants.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <sstream>
#include <utility>

namespace labw::art_modern {

void ArtRead::generate_pairwise_aln()
{
    hts_pos_t pos_on_aln_str = 0;
    hts_pos_t pos_on_query = 0;
    hts_pos_t pos_on_ref = 0;
    const std::size_t maxk = art_params_.read_len + 1 + 1 + indel_.size();
    aln_query_.resize(maxk);
    aln_ref_.resize(maxk);
#if (1)
    hts_pos_t num_match = 0;
    for (const auto& [this_pos_on_aln_str, indel] : indel_) {
        num_match = this_pos_on_aln_str - pos_on_aln_str;
        std::memcpy(aln_query_.data() + pos_on_aln_str, query_.data() + pos_on_query, num_match);
        std::memcpy(aln_ref_.data() + pos_on_aln_str, ref_.data() + pos_on_ref, num_match);
        pos_on_query += num_match;
        pos_on_ref += num_match;
        pos_on_aln_str = this_pos_on_aln_str;
        if (indel == ALN_GAP) {
            aln_query_[pos_on_aln_str] = ALN_GAP;
            aln_ref_[pos_on_aln_str] = ref_[pos_on_ref];
            pos_on_ref++;
        } else {
            aln_query_[pos_on_aln_str] = query_[pos_on_query];
            aln_ref_[pos_on_aln_str] = ALN_GAP;
            pos_on_query++;
        }
        pos_on_aln_str++;
    }
    num_match = static_cast<hts_pos_t>(ref_.size() - pos_on_ref);
#ifdef CEU_CM_IS_DEBUG
    if (num_match < 0) {
        BOOST_LOG_TRIVIAL(fatal) << "num_match < 0: " << num_match;
        abort_mpi();
    }
#endif
    std::memcpy(aln_query_.data() + pos_on_aln_str, query_.data() + pos_on_query, num_match);
    std::memcpy(aln_ref_.data() + pos_on_aln_str, ref_.data() + pos_on_ref, num_match);
    pos_on_aln_str += num_match;
#else // Old version for historical purposes
    while (pos_on_ref < seq_ref_.size()) {
        const auto find = indel_.find(pos_on_aln_str);
        if (find == indel_.end()) { // No indel
            aln_read_[pos_on_aln_str] = seq_read_[pos_on_read];
            aln_ref_[pos_on_aln_str] = seq_ref_[pos_on_ref];
            pos_on_read++;
            pos_on_ref++;
        } else if (find->second == ALN_GAP) { // Deletion
            aln_read_[pos_on_aln_str] = ALN_GAP;
            aln_ref_[pos_on_aln_str] = seq_ref_[pos_on_ref];
            pos_on_ref++;
        } else { // Insertion
            aln_read_[pos_on_aln_str] = find->second;
            aln_ref_[pos_on_aln_str] = ALN_GAP;
            pos_on_read++;
        }
        pos_on_aln_str++;
    }
    while (true) { // Insertions after reference
        const auto find = indel_.find(pos_on_aln_str);
        if (find == indel_.end()) {
            break;
        }
        aln_read_[pos_on_aln_str] = find->second;
        aln_ref_[pos_on_aln_str] = ALN_GAP;
        pos_on_aln_str++;
    }
#endif
    aln_query_.resize(pos_on_aln_str);
    aln_ref_.resize(pos_on_aln_str);

#ifdef CEU_CM_IS_DEBUG
    if (aln_query_.size() != aln_ref_.size()) {
        BOOST_LOG_TRIVIAL(fatal) << "Aligned read size (" << aln_query_.size() << ") != aligned ref size ("
                                 << aln_ref_.size() << ")";
        except_();
    }
    auto reconst_ref = aln_ref_;
    auto reconst_query = aln_query_;
    reconst_ref.erase(std::remove(reconst_ref.begin(), reconst_ref.end(), ALN_GAP), reconst_ref.end());
    reconst_query.erase(std::remove(reconst_query.begin(), reconst_query.end(), ALN_GAP), reconst_query.end());

    if (reconst_ref != ref_) {
        BOOST_LOG_TRIVIAL(error) << "Reconstructed reference != reference:";
        except_();
    }
    if (reconst_query != query_) {
        BOOST_LOG_TRIVIAL(error) << "Reconstructed query != query:";
        except_();
    }
#endif
}

void ArtRead::generate_snv_on_qual(const bool is_first_read)
{
    if (!art_params_.sep_flag) {
        art_params_.qdist.get_read_qual(qual_, rprob_, is_first_read);
    } else if (is_first_read) {
        art_params_.qdist.get_read_qual_sep_1(qual_, query_, rprob_);
    } else {
        art_params_.qdist.get_read_qual_sep_2(qual_, query_, rprob_);
    }
    char achar = 0;
    rprob_.r_probs();
    for (decltype(qual_.size()) i = 0; i < qual_.size(); i++) {
        if (query_[i] == 'N') {
            qual_[i] = MIN_QUAL;
            continue;
        }
        if (rprob_.tmp_probs_[i] < art_params_.err_prob[qual_[i]]) {
            do {
                achar = rprob_.rand_base();
            } while (query_[i] == achar);
            query_[i] = achar;
        }
    }
}

int ArtRead::generate_indels(const bool is_read_1)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    int i = 0;
    int j = 0;
    int pos = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
    // deletion
    rprob_.r_probs(per_base_del_rate.size());
    for (i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= rprob_.tmp_probs_[i]) {
            del_len = i + 1;
            j = i;
            while (j >= 0) {
                pos = rprob_.rand_pos_on_read_not_head_and_tail();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }
    rprob_.r_probs(per_base_ins_rate.size());
    for (i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if (art_params_.read_len - del_len - ins_len < i + 1) {
            continue; // ensure that there are stil enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= rprob_.tmp_probs_[i]) {
            ins_len = i + 1;
            j = i;
            while (j >= 0) {
                pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }
    return ins_len - del_len;
}

int ArtRead::generate_indels_2(const bool is_read_1)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    int pos = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

    rprob_.r_probs(per_base_ins_rate.size());
    for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= rprob_.tmp_probs_[i]) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }

    // deletion
    rprob_.r_probs(per_base_del_rate.size());
    for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if (art_params_.read_len - del_len - ins_len < i + 1) {
            continue; // ensure that enough unchanged position for mutation
        }

        if (per_base_del_rate[i] >= rprob_.tmp_probs_[i]) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                pos = rprob_.rand_pos_on_read_not_head_and_tail();
                if (pos == 0) {
                    continue;
                }
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = ALN_GAP;
                    j--;
                }
            }
            break;
        }
    }
    return ins_len - del_len;
}

void ArtRead::ref2read(std::string seq_ref, const bool is_plus_strand, const hts_pos_t pos_on_contig)
{
    pos_on_contig_ = pos_on_contig;
    is_plus_strand_ = is_plus_strand;
    ref_ = std::move(seq_ref);
    normalize_inplace(ref_);
    if (!is_plus_strand) {
        revcomp_inplace(ref_);
    }
    hts_pos_t pos_on_aln_str = 0;
    hts_pos_t pos_on_read = 0;
    hts_pos_t pos_on_ref = 0;
    const std::size_t maxk = art_params_.read_len + 1 + 1 + indel_.size();
    aln_query_.resize(maxk);
    aln_ref_.resize(maxk);
#if (1)
    hts_pos_t num_match = 0;

    for (const auto& [this_pos_on_aln_str, indel] : indel_) {
        num_match = this_pos_on_aln_str - pos_on_aln_str;
        std::memcpy(query_.data() + pos_on_read, ref_.data() + pos_on_ref, num_match);
        pos_on_read += num_match;
        pos_on_ref += num_match;
        pos_on_aln_str = this_pos_on_aln_str;
        if (indel == ALN_GAP) {
            pos_on_ref++;
        } else {
            query_[pos_on_read] = indel;
            pos_on_read++;
        }
        pos_on_aln_str++;
    }
    num_match = static_cast<hts_pos_t>(ref_.size() - pos_on_ref);

#ifdef CEU_CM_IS_DEBUG
    if (num_match < 0) {
        BOOST_LOG_TRIVIAL(fatal) << "num_match < 0: " << num_match;
        abort_mpi();
    }
#endif
    std::memcpy(query_.data() + pos_on_read, ref_.data() + pos_on_ref, num_match);
#else // Old code for historical purposes
    for (decltype(seq_ref_.size()) pos_on_ref = 0; pos_on_ref < seq_ref_.size();) {
        const auto find = indel_.find(k);
        if (find == indel_.end()) { // No indel
            seq_read_[pos_on_read] = seq_ref_[pos_on_ref];
            pos_on_ref++;
            pos_on_read++;
        } else if (find->second == ALN_GAP) { // Deletion
            pos_on_ref++;
        } else { // Insertion
            seq_read_[pos_on_read] = find->second;
            pos_on_read++;
        }
        k++;
    }
    while (true) { // Insertions after reference
        const auto find = indel_.find(k);
        if (find == indel_.end()) {
            break;
        }
        seq_read_[pos_on_read] = find->second;
        pos_on_read++;
        k++;
    }
#endif
#ifdef CEU_CM_IS_DEBUG
    if (static_cast<int>(query_.size()) != art_params_.read_len) {
        BOOST_LOG_TRIVIAL(error) << "Generated read length (" << query_.size()
                                 << ") is not equal to designed read length (" << art_params_.read_len << ")";
        except_();
    }
#endif
}

void ArtRead::except_() const
{
    BOOST_LOG_TRIVIAL(error) << "Read   : " << read_name_;
    BOOST_LOG_TRIVIAL(error) << "Contig : " << contig_name_ << ":" << pos_on_contig_ << ":"
                             << (is_plus_strand_ ? "+" : "-");
    BOOST_LOG_TRIVIAL(error) << "Query  : " << query_;
    BOOST_LOG_TRIVIAL(error) << "Ref    : " << ref_;
    BOOST_LOG_TRIVIAL(error) << "Qual   : " << qual_to_str(qual_);
    BOOST_LOG_TRIVIAL(error) << "AQuery : " << aln_query_;
    BOOST_LOG_TRIVIAL(error) << "ARef   : " << aln_ref_;
    std::ostringstream indel_oss;
    indel_oss << "{";
    for (const auto& [pos, indel] : indel_) {
        indel_oss << pos << ":" << indel << ", ";
    }
    indel_oss << "}";
    BOOST_LOG_TRIVIAL(error) << "Indel  : " << indel_oss.str();
    abort_mpi();
}

ArtRead::ArtRead(const ArtParams& art_params, std::string contig_name, std::string read_name, Rprob& rprob)
    : art_params_(art_params)
    , contig_name_(std::move(contig_name))
    , read_name_(std::move(read_name))
    , rprob_(rprob)
{
    query_.resize(art_params_.read_len);
    qual_.resize(art_params_.read_len);
}

PairwiseAlignment ArtRead::to_pwa()
{
    auto qual_str = qual_to_str(qual_);
    return { std::move(read_name_), std::move(contig_name_), std::move(query_), std::move(ref_), std::move(qual_str),std::move(qual_),
        std::move(aln_query_), std::move(aln_ref_), pos_on_contig_, is_plus_strand_ };
}
bool ArtRead::is_good() const { return std::count(query_.begin(), query_.end(), 'N') <= art_params_.max_n; }

} // namespace labw::art_modern
