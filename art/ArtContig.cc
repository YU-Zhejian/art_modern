#include "ArtContig.hh"
#include "art_modern_constants.hh"
#include "random_generator.hh"
namespace labw::art_modern {

/**
 * SE:@code
                |----------->
             ------------------------------------
 OR
             ------------------------------------
                               <-----------|
 * @endcode
 * @param is_plus_strand
 * @param read_1
 */
void ArtContig::generate_read_se(const bool is_plus_strand, ArtRead& read_1)
{
    read_1.is_plus_strand = is_plus_strand;
    const auto pos_1 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        ? 0
        : rprob_.randint(0, static_cast<int>(valid_region_) + 1);
    auto slen_1 = read_1.generate_indels(true);
    // ensure get a fixed read length
    if (pos_1 + art_params_.read_len - slen_1 > ref_len_) {
        slen_1 = read_1.generate_indels_2(true);
    }
    const auto ref = normalize(fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1));
    read_1.seq_ref = read_1.is_plus_strand ? ref : revcomp(ref);
    read_1.pos_on_contig = pos_1;
    read_1.ref2read();
}

/**
 *
 * PE: @code
                |----------->
             ------------------------------------
                               <-----------|

 * @endcode
 * MP: @code
                <-----------|
             ------------------------------------
                               |----------->

 * @endcode
 *
 * @param is_plus_strand
 * @param is_mp
 * @param read_1
 * @param read_2
 */
void ArtContig::generate_read_pe(const bool is_plus_strand, const bool is_mp, ArtRead& read_1, ArtRead& read_2)
{
    const hts_pos_t fragment_len = generate_fragment_length();
    const hts_pos_t fragment_start = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        ? 0
        : rprob_.randint(0, static_cast<int>(ref_len_ - fragment_len) + 1);
    const hts_pos_t fragment_end = fragment_start + fragment_len;

    const hts_pos_t pos_1 = is_mp == is_plus_strand ? fragment_end - art_params_.read_len : fragment_start;
    const hts_pos_t pos_2 = is_mp == is_plus_strand ? fragment_start : fragment_end - art_params_.read_len;

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
    const auto ref1 = normalize(fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1));
    const auto ref2 = normalize(fasta_fetch_->fetch(seq_id_, pos_2, pos_2 + art_params_.read_len - slen_2));

    if (is_plus_strand) {
        read_1.seq_ref = ref1;
        read_2.seq_ref = revcomp(ref2);
    } else {
        read_1.seq_ref = revcomp(ref1);
        read_2.seq_ref = ref2;
    }

    read_1.pos_on_contig = pos_1;
    read_1.ref2read();
    read_2.pos_on_contig = pos_2;
    read_2.ref2read();
}

ArtContig::ArtContig(BaseFastaFetch* fasta_fetch, const size_t seq_id, const ArtParams& art_params, Rprob& rprob)
    : seq_name(fasta_fetch->seq_name(seq_id))
    , fasta_fetch_(fasta_fetch)
    , art_params_(art_params)
    , rprob_(rprob)
    , seq_id_(seq_id)
    , valid_region_(fasta_fetch->seq_len(seq_id) - art_params.read_len)
    , ref_len_(fasta_fetch->seq_len(seq_id))
{
}
hts_pos_t ArtContig::generate_fragment_length() const
{
    hts_pos_t fragment_len;
    if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        || art_params_.pe_dist_mean_minus_2_std > ref_len_) {
        // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to
        // be reference length
        fragment_len = ref_len_;
    } else {
        fragment_len = 0;
        while (fragment_len < art_params_.read_len || fragment_len > ref_len_) {
            fragment_len = rprob_.insertion_length();
        }
    }
    return fragment_len;
}

}