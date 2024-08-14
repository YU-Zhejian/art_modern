#include <string>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <utility>

#include "ArtContig.hh"
#include "Empdist.hh"

using namespace std;
using namespace labw::art_modern;

ArtRead ArtContig::generate_read_se() const
{
    ArtRead read_1(_art_params);
    auto ref_len = static_cast<int>(_ref_seq.length());
    auto pos_1 = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE
        ? 0
        : static_cast<int>(floor(rprob.r_prob() * _valid_region));
    auto slen_1 = read_1.generate_indels(_art_params.read_len, true);
    // ensure get a fixed read length
    if (pos_1 + _art_params.read_len - slen_1 > ref_len) {
        slen_1 = read_1.generate_indels_2(_art_params.read_len, true);
    }
    read_1.is_plus_strand = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE || rprob.r_prob() < 0.5;
    if (read_1.is_plus_strand) {
        //    |----------->
        // ------------------------------------
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);
    } else {
        // ------------------------------------
        //                   <-----------|
        read_1.seq_ref = revcomp(_ref_seq.substr(_valid_region - pos_1, _art_params.read_len - slen_1));
    }
    read_1.bpos = pos_1;
    read_1.ref2read();
    return read_1;
}

// matepair-end read: the second read is reverse complemenaty strand
ArtReadPair ArtContig::generate_read_mp() const
{
    ArtRead read_1(_art_params);
    ArtRead read_2(_art_params);
    int fragment_len;
    auto ref_len = static_cast<int>(_ref_seq.length());
    if (_art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE) {
        fragment_len = ref_len;
    } else {
        auto rng = boost::mt19937();
        auto gaussian = boost::normal_distribution<>(_art_params.pe_frag_dist_mean, _art_params.pe_frag_dist_std_dev);
        if (_art_params.pe_frag_dist_mean - 2 * _art_params.pe_frag_dist_std_dev > ref_len) {
            // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to
            // be reference length
            fragment_len = ref_len;
        } else {
            fragment_len = 0;
            while (fragment_len < _art_params.read_len || fragment_len > _ref_seq.length()) {
                fragment_len = static_cast<int>(gaussian(rng));
            }
        }
    }

    auto pos_1 = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE
        ? ref_len - _art_params.read_len
        : static_cast<int>(floor((ref_len - fragment_len) * rprob.r_prob()) + fragment_len - _art_params.read_len);
    auto pos_2 = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE
        ? ref_len - _art_params.read_len
        : ref_len - (pos_1 + 2 * _art_params.read_len - fragment_len);
    // Exceptions at this step for unknown reasons.
    auto slen_1 = read_1.generate_indels(_art_params.read_len, true);
    auto slen_2 = read_2.generate_indels(_art_params.read_len, false);

    // ensure get a fixed read length
    if ((pos_1 + _art_params.read_len - slen_1) > ref_len) {
        slen_1 = read_1.generate_indels_2(_art_params.read_len, true);
    }
    if ((pos_2 + _art_params.read_len - slen_2) > ref_len) {
        slen_2 = read_2.generate_indels_2(_art_params.read_len, false);
    }

    bool is_plus_strand = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE || rprob.r_prob() < 0.5;
    if (is_plus_strand) {
        // R1                  |----------->
        //   ------------------------------------
        // R2   <-----------|
        read_1.is_plus_strand = true;
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);

        read_2.is_plus_strand = false;
        read_2.seq_ref = revcomp(_ref_seq.substr(_valid_region - pos_2, _art_params.read_len - slen_2));
    } else {
        // R2   <-----------|
        //   ------------------------------------
        // R1                  |----------->
        read_1.is_plus_strand = false;
        read_1.seq_ref = revcomp(_ref_seq.substr(_valid_region - pos_1, _art_params.read_len - slen_1));
        read_2.is_plus_strand = true;
        read_2.seq_ref = _ref_seq.substr(pos_2, _art_params.read_len - slen_2);
    }
    read_1.bpos = pos_1;
    read_1.ref2read();
    read_2.bpos = pos_2;
    read_2.ref2read();
    return { read_1, read_2 };
}

// paired-end read: the second read is reverse complemenaty strand
ArtReadPair ArtContig::generate_read_pe() const
{
    ArtRead read_1(_art_params);
    ArtRead read_2(_art_params);
    int fragment_len;
    auto ref_len = static_cast<int>(_ref_seq.length());
    if (_art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE) {
        fragment_len = ref_len;
    } else {
        auto rng = boost::mt19937();
        // FIXME: _art_params.pe_frag_dist_std_dev seems not working
        auto gaussian = boost::normal_distribution<>(_art_params.pe_frag_dist_mean, _art_params.pe_frag_dist_std_dev);
        if (_art_params.pe_frag_dist_mean - 2 * _art_params.pe_frag_dist_std_dev > ref_len) {
            // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to
            // be reference length
            fragment_len = ref_len;
        } else {
            fragment_len = 0;
            while (fragment_len < _art_params.read_len || fragment_len > _ref_seq.length()) {
                fragment_len = static_cast<int>(gaussian(rng));
            }
        }
    }
    auto pos_1 = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE
        ? 0
        : static_cast<int>(floor((ref_len - fragment_len) * rprob.r_prob()));
    auto pos_2 = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE ? 0 : ref_len - pos_1 - fragment_len;
    int slen_1 = read_1.generate_indels(_art_params.read_len, true);
    int slen_2 = read_2.generate_indels(_art_params.read_len, false);

    // ensure get a fixed read length
    if ((pos_1 + _art_params.read_len - slen_1) > ref_len) {
        slen_1 = read_1.generate_indels_2(_art_params.read_len, true);
    }
    if ((pos_2 + _art_params.read_len - slen_2) > ref_len) {
        slen_2 = read_2.generate_indels_2(_art_params.read_len, false);
    }

    bool is_plus_strand = _art_params.art_simulation_mode == ART_SIMULATION_MODE::TEMPLATE || rprob.r_prob() < 0.5;

    if (is_plus_strand) {
        //    |----------->
        // ------------------------------------
        //                   <-----------|
        read_1.is_plus_strand = true;
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);
        read_2.is_plus_strand = false;
        read_2.seq_ref = revcomp(_ref_seq.substr(_valid_region - pos_2, _art_params.read_len - slen_2));
    } else {
        //                   <-----------|
        // ------------------------------------
        //    |----------->
        read_1.is_plus_strand = false;
        read_1.seq_ref = revcomp(_ref_seq.substr(_valid_region - pos_1, _art_params.read_len - slen_1));
        read_2.is_plus_strand = true;
        read_2.seq_ref = _ref_seq.substr(pos_2, _art_params.read_len - slen_2);
    }
    read_1.bpos = pos_1;
    read_1.ref2read();
    read_2.bpos = pos_2;
    read_2.ref2read();
    return { read_1, read_2 };
}

ArtContig::ArtContig(std::string ref_seq, std::string id, const ArtParams& art_params)
    : _art_params(art_params)
    , _ref_seq(ref_seq)
    , _id(std::move(id))
{
    std::replace(ref_seq.begin(), ref_seq.end(), 'U', 'T');

    const auto ref_len = static_cast<int>(ref_seq.length());
    boost::algorithm::to_upper(ref_seq);
    _valid_region = ref_len - _art_params.read_len;
}
