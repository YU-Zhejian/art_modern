#include "ArtRead.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"

using namespace std;

namespace labw::art_modern {

void ArtRead::generate_pairwise_aln()
{
    int k = 0;
    std::size_t pos_on_read = 0;
    std::size_t maxk = art_params_.read_len + 1 + 1 + indel_.size();
    aln_read_.resize(maxk);
    aln_ref_.resize(maxk);

    for (decltype(seq_ref_.size()) pos_on_ref = 0; pos_on_ref < seq_ref_.size();) {
        if (indel_.find(k) == indel_.end()) { // No indel
            aln_read_[k] = seq_read_[pos_on_read];
            aln_ref_[k] = seq_ref_[pos_on_ref];
            pos_on_read++;
            pos_on_ref++;
        } else if (indel_[k] == ALN_GAP) { // Deletion
            aln_read_[k] = ALN_GAP;
            aln_ref_[k] = seq_ref_[pos_on_ref];
            pos_on_ref++;
        } else { // Insertion
            aln_read_[k] = indel_[k];
            aln_ref_[k] = ALN_GAP;
            pos_on_read++;
        }
        k++;
    }
    while (indel_.find(k) != indel_.end()) { // Insertions after reference
        aln_read_[k] = indel_[k];
        aln_ref_[k] = ALN_GAP;
        pos_on_read++;
        k++;
    }
    aln_read_.resize(k);
    aln_ref_.resize(k);
}

void ArtRead::generate_snv_on_qual(const bool is_first_read, std::vector<double>& tmp_qual_probs)
{
    if (!art_params_.sep_flag) {
        art_params_.qdist.get_read_qual(qual_, art_params_.read_len, rprob_, is_first_read);
    } else if (is_first_read) {
        art_params_.qdist.get_read_qual_sep_1(qual_, seq_read_, rprob_);
    } else {
        art_params_.qdist.get_read_qual_sep_2(qual_, seq_read_, rprob_);
    }
    char achar;
    rprob_.r_probs(tmp_qual_probs);
    for (decltype(qual_.size()) i = 0; i < qual_.size(); i++) {
        if (seq_read_[i] == 'N') {
            qual_[i] = MIN_QUAL;
            continue;
        }
        if (tmp_qual_probs[i] < art_params_.err_prob[qual_[i]]) {
            achar = seq_read_[i];
            while (seq_read_[i] == achar) {
                achar = rprob_.rand_base();
            }
            seq_read_[i] = achar;
        }
    }
}

int ArtRead::generate_indels(const bool is_read_1, std::vector<double>& probs_indel)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    int i;
    int j;
    int pos;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;
    // deletion
    probs_indel.resize(per_base_del_rate.size());
    rprob_.r_probs(probs_indel);
    for (i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (per_base_del_rate[i] >= probs_indel[i]) {
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
    probs_indel.resize(per_base_ins_rate.size());
    rprob_.r_probs(probs_indel);
    for (i = static_cast<int>(per_base_ins_rate.size() - 1); i >= 0; i--) {
        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }
        if (per_base_ins_rate[i] >= probs_indel[i]) {
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
    return (ins_len - del_len);
};

int ArtRead::generate_indels_2(const bool is_read_1, std::vector<double>& probs_indel)
{
    indel_.clear();
    int ins_len = 0;
    int del_len = 0;
    const auto& per_base_del_rate = is_read_1 ? art_params_.per_base_del_rate_1 : art_params_.per_base_del_rate_2;
    const auto& per_base_ins_rate = is_read_1 ? art_params_.per_base_ins_rate_1 : art_params_.per_base_ins_rate_2;

    probs_indel.resize(per_base_ins_rate.size());
    rprob_.r_probs(probs_indel);
    for (auto i = static_cast<int>(per_base_ins_rate.size()) - 1; i >= 0; i--) {
        if (per_base_ins_rate[i] >= probs_indel[i]) {
            ins_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = rprob_.rand_pos_on_read();
                if (indel_.find(pos) == indel_.end()) {
                    indel_[pos] = rprob_.rand_base();
                    j--;
                }
            }
            break;
        }
    }

    // deletion
    probs_indel.resize(per_base_ins_rate.size());
    rprob_.r_probs(probs_indel);
    for (auto i = static_cast<int>(per_base_del_rate.size()) - 1; i >= 0; i--) {
        if (del_len == ins_len) {
            break;
        }

        if ((art_params_.read_len - del_len - ins_len) < (i + 1)) {
            continue; // ensure that enough unchanged position for mutation
        }

        if (per_base_del_rate[i] >= probs_indel[i]) {
            del_len = i + 1;
            for (int j = i; j >= 0;) {
                auto pos = rprob_.rand_pos_on_read_not_head_and_tail();
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
    return (ins_len - del_len);
};

void ArtRead::ref2read(std::string seq_ref, const bool is_plus_strand, const hts_pos_t pos_on_contig)
{
    pos_on_contig_ = pos_on_contig;
    is_plus_strand_ = is_plus_strand;
    seq_ref_ = std::move(seq_ref);
    normalize_inplace(seq_ref_);
    if (!is_plus_strand) {
        revcomp_inplace(seq_ref_);
    }

    int k = 0;
    int pos_on_read = 0;
    for (decltype(seq_ref_.size()) pos_on_ref = 0; pos_on_ref < seq_ref_.size();) {
        if (indel_.find(k) == indel_.end()) { // No indel
            seq_read_[pos_on_read] = seq_ref_[pos_on_ref];
            pos_on_ref++;
            pos_on_read++;
        } else if (indel_[k] == ALN_GAP) { // Deletion
            pos_on_ref++;
        } else { // Insertion
            seq_read_[pos_on_read] = indel_[k];
            pos_on_read++;
        }
        k++;
    }
    while (indel_.find(k) != indel_.end()) { // Insertions after reference
        seq_read_[pos_on_read] = indel_[k];
        pos_on_read++;
        k++;
    }
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        throw ReadGenerationException("Generated seq_read_ with unequal sizes");
    }
}

ArtRead::ArtRead(
    const ArtParams& art_params, const std::string& contig_name, const std::string& read_name, Rprob& rprob)
    : art_params_(art_params)
    , contig_name_(contig_name)
    , read_name_(read_name)
    , rprob_(rprob)
{
    seq_read_.resize(art_params_.read_len);
    qual_.resize(art_params_.read_len);
}

PairwiseAlignment ArtRead::to_pwa()
{
    std::string qual_str = qual_to_str(qual_);
    return { std::move(read_name_), std::move(contig_name_), std::move(seq_read_), std::move(seq_ref_), qual_str,
        std::move(aln_read_), std::move(aln_ref_), pos_on_contig_, is_plus_strand_ };
}
bool ArtRead::is_good() const
{
    if (std::count(seq_read_.begin(), seq_read_.end(), 'N') > 0) { // TODO: Add params back.
        return false;
    }
    if (static_cast<int>(seq_read_.size()) != art_params_.read_len) {
        goto error;
    }
    if (static_cast<int>(qual_.size()) != art_params_.read_len) {
        goto error;
    }
    if (aln_read_.size() != aln_ref_.size()) {
        goto error;
    }
    return true;
error:
    // #ifdef CEU_CM_IS_DEBUG
    //         abort_mpi();
    // #else
    return false; // FIXME: No idea why this occurs.
    // #endif
}

} // namespace labw::art_modern; // namespace labw