#include <string>

#include "ArtContig.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"
namespace labw::art_modern {

void ArtContig::generate_read_se(const bool is_plus_strand, ArtRead& read_1)
{
    read_1.is_plus_strand = is_plus_strand;

    const auto pos_1 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        ? 0
        : static_cast<hts_pos_t>(std::floor(rprob_.r_prob() * static_cast<double>(valid_region_)));
    auto slen_1 = read_1.generate_indels(true);
    // ensure get a fixed read length
    if (pos_1 + art_params_.read_len - slen_1 > ref_len_) {
        slen_1 = read_1.generate_indels_2(true);
    }
    const auto ref = normalize(fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1));
    read_1.seq_ref = read_1.is_plus_strand ? ref : revcomp(ref);
    read_1.bpos = pos_1;
    read_1.ref2read();
    if (read_1.seq_read.size() != art_params_.read_len) {
        throw ReadGenerationException("Generated seq_read with unequal sizes");
    }
}

// paired-end read: the second read is reverse complemenaty strand
void ArtContig::generate_read_pe(const bool is_plus_strand, const bool is_mp, ArtRead& read_1, ArtRead& read_2)
{
    hts_pos_t fragment_len;
    if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
        fragment_len = ref_len_;
    } else {
        if (art_params_.pe_dist_mean_minus_2_std > ref_len_) {
            // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to
            // be reference length
            fragment_len = ref_len_;
        } else {
            fragment_len = 0;
            while (fragment_len < art_params_.read_len || fragment_len > ref_len_) {
                fragment_len = rprob_.insertion_length();
            }
        }
    }
    hts_pos_t pos_1;
    hts_pos_t pos_2;

    if (is_mp) {
        pos_1 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? ref_len_ - art_params_.read_len
            : static_cast<hts_pos_t>(std::floor(static_cast<double>(ref_len_ - fragment_len) * rprob_.r_prob()))
                + fragment_len - art_params_.read_len;
        pos_2 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? ref_len_ - art_params_.read_len
            : ref_len_ - (pos_1 + 2 * art_params_.read_len - fragment_len);
    } else {
        pos_1 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? 0
            : static_cast<hts_pos_t>(std::floor(static_cast<double>(ref_len_ - fragment_len) * rprob_.r_prob()));
        pos_2 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE ? 0 : ref_len_ - pos_1 - fragment_len;
    }
    int slen_1 = read_1.generate_indels(true);
    int slen_2 = read_2.generate_indels(false);

    // ensure get a fixed read length
    if ((pos_1 + art_params_.read_len - slen_1) > ref_len_) {
        slen_1 = read_1.generate_indels_2(true);
    }
    if ((pos_2 + art_params_.read_len - slen_2) > ref_len_) {
        slen_2 = read_2.generate_indels_2(false);
    }
    read_1.is_plus_strand = is_plus_strand;
    read_2.is_plus_strand = !is_plus_strand;

    if (is_plus_strand) {
        read_1.seq_ref = normalize(fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1));
        read_2.seq_ref = revcomp(normalize(fasta_fetch_->fetch(
            seq_id_, valid_region_ - pos_2, valid_region_ - pos_2 + art_params_.read_len - slen_2)));
    } else {
        read_1.seq_ref = revcomp(normalize(fasta_fetch_->fetch(
            seq_id_, valid_region_ - pos_1, valid_region_ - pos_1 + art_params_.read_len - slen_1)));
        read_2.seq_ref = normalize(fasta_fetch_->fetch(seq_id_, pos_2, pos_2 + art_params_.read_len - slen_2));
    }
    read_1.bpos = pos_1;
    read_1.ref2read();
    read_2.bpos = pos_2;
    read_2.ref2read();
    if (read_1.seq_read.size() != art_params_.read_len && read_2.seq_read.size() != art_params_.read_len) {
        throw ReadGenerationException("Generated seq_read with unequal sizes");
    }
}

ArtContig::ArtContig(BaseFastaFetch* fasta_fetch, const size_t seq_id, const ArtParams& art_params, Rprob& rprob)
    : seq_name_(fasta_fetch->seq_name(seq_id))
    , ref_len_(fasta_fetch->seq_len(seq_id))
    , fasta_fetch_(fasta_fetch)
    , art_params_(art_params)
    , rprob_(rprob)
    , seq_id_(seq_id)
    , valid_region_(fasta_fetch->seq_len(seq_id) - art_params.read_len)
{
}

}