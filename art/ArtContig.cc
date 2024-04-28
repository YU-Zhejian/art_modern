#include <string>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <utility>

#include "ArtContig.hh"
#include "Empdist.hh"

using namespace std;
using namespace labw::art_modern;

void ArtContig::mask_n_region(int max_num_n)
{
    _masked_pos.clear(); // reset mask position
    if (_ref_seq.size() < _art_params.read_len)
        return;
    int num_n = 0;
    for (int i = 0; i < _art_params.read_len - 1; i++) {
        if (_ref_seq[i] == 'N')
            num_n++;
    }
    for (size_t i = 0; i <= _ref_seq.size() - _art_params.read_len; i++) {
        if (_ref_seq[i + _art_params.read_len - 1] == 'N')
            num_n++;
        if (num_n >= max_num_n) {
            _masked_pos.insert(i);
        }
        if (_ref_seq[i] == 'N')
            num_n--;
    }
}

ArtRead ArtContig::generate_read_se(bool is_amplicon) const
{
    ArtRead read_1(_art_params);
    auto ref_len = static_cast<int>(_ref_seq.length());
    auto pos_1 = is_amplicon ? 0 : static_cast<int>(floor(r_prob() * _valid_region));
    auto slen_1 = read_1.generate_indels(_art_params.read_len, true);
    // ensure get a fixed read length
    if (pos_1 + _art_params.read_len - slen_1 > ref_len) {
        slen_1 = read_1.generate_indels_2(_art_params.read_len, true);
    }
    read_1.is_plus_strand = is_amplicon || r_prob() < 0.5;
    // string read_info_1;
    // int aln_start;
    // int aln_end;
    if (read_1.is_plus_strand) {
        //    |----------->
        // ------------------------------------
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % pos_1 % (_art_params.read_len - slen_1 + pos_1) % "+").str();
        // aln_start = pos_1;
        // aln_end = _art_params.read_len - slen_1 + pos_1;
    } else {
        // ------------------------------------
        //                   <-----------|
        read_1.seq_ref = _ref_seq_cmp.substr(pos_1, _art_params.read_len - slen_1);
        // FIXME: Whether those plus ones are needed are not sure.
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % (ref_len - pos_1 + 1) % (ref_len - pos_1 + 1 - _art_params.read_len + slen_1) % "-").str();
        // aln_start = ref_len - pos_1 + 1 - _art_params.read_len + slen_1;
        // aln_end = ref_len - pos_1 + 1;
    }
    // std::cout << "R1: " << read_info_1 << " Dist: " << aln_end - aln_start + 1 << endl;
    read_1.bpos = pos_1;
    read_1.ref2read();
    return read_1;
}

// matepair-end read: the second read is reverse complemenaty strand
ArtReadPair ArtContig::generate_read_mp(bool is_amplicon) const
{
    ArtRead read_1(_art_params);
    ArtRead read_2(_art_params);
    int fragment_len;
    auto ref_len = static_cast<int>(_ref_seq.length());
    if (is_amplicon) {
        fragment_len = ref_len;
    } else {
        auto rng = boost::mt19937();
        auto gaussian = boost::normal_distribution<>(_art_params.pe_frag_dist_mean, _art_params.pe_frag_dist_std_dev);
        if (_art_params.pe_frag_dist_mean - 2 * _art_params.pe_frag_dist_std_dev > ref_len) {
            // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to be reference length
            fragment_len = ref_len;
        } else {
            fragment_len = 0;
            while (fragment_len < _art_params.read_len || fragment_len > _ref_seq.length()) {
                fragment_len = static_cast<int>(gaussian(rng));
            }
        }
    }

    auto pos_1 = is_amplicon ? ref_len - _art_params.read_len : static_cast<int>(floor((ref_len - fragment_len) * r_prob()) + fragment_len - _art_params.read_len);
    auto pos_2 = is_amplicon ? ref_len - _art_params.read_len : ref_len - (pos_1 + 2 * _art_params.read_len - fragment_len);
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

    bool is_plus_strand = is_amplicon || r_prob() < 0.5;
    // string read_info_1;
    // string read_info_2;
    // int aln_end;
    // int aln_start;
    if (is_plus_strand) {
        // R1                  |----------->
        //   ------------------------------------
        // R2   <-----------|
        read_1.is_plus_strand = true;
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % pos_1 % (_art_params.read_len - slen_1 + pos_1) % "+").str();
        // aln_end = pos_1 + _art_params.read_len - slen_1;

        read_2.is_plus_strand = false;
        read_2.seq_ref = _ref_seq_cmp.substr(pos_2, _art_params.read_len - slen_2);
        // read_info_2 = (boost::format("%s:%d-%d:%s") % _id % (ref_len - pos_2) % (ref_len - pos_2 - _art_params.read_len + slen_2) % "-").str();
        // aln_start = ref_len - pos_2 - _art_params.read_len + slen_2;
    } else {
        // R2   <-----------|
        //   ------------------------------------
        // R1                  |----------->
        read_1.is_plus_strand = false;
        read_1.seq_ref = _ref_seq_cmp.substr(pos_1, _art_params.read_len - slen_1);
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % pos_1 % (pos_1 + _art_params.read_len - slen_1) % "-").str();
        // aln_end = pos_1 + _art_params.read_len - slen_1;

        read_2.is_plus_strand = true;
        read_2.seq_ref = _ref_seq.substr(pos_2, _art_params.read_len - slen_2);
        // read_info_2 = (boost::format("%s:%d-%d:%s") % _id % (ref_len - pos_2) % (ref_len - pos_2 - _art_params.read_len + slen_2) % "+").str();
        // aln_start = ref_len - pos_2 - _art_params.read_len + slen_2;
    }
    // std::cout << "R1: " << read_info_1 << ", R2: " << read_info_2 << " Dist: " << aln_end - aln_start + 1 << endl;
    read_1.bpos = pos_1;
    read_1.ref2read();
    read_2.bpos = pos_2;
    read_2.ref2read();
    // On reverse: pos -= del_len - ins_len; // FIXME: adjust start position
    return { read_1, read_2 };
}

// paired-end read: the second read is reverse complemenaty strand
ArtReadPair ArtContig::generate_read_pe(bool is_amplicon) const
{
    ArtRead read_1(_art_params);
    ArtRead read_2(_art_params);
    int fragment_len;
    auto ref_len = static_cast<int>(_ref_seq.length());
    if (is_amplicon) {
        fragment_len = ref_len;
    } else {
        auto rng = boost::mt19937();
        // FIXME: _art_params.pe_frag_dist_std_dev seems not working
        auto gaussian = boost::normal_distribution<>(_art_params.pe_frag_dist_mean, _art_params.pe_frag_dist_std_dev);
        if (_art_params.pe_frag_dist_mean - 2 * _art_params.pe_frag_dist_std_dev > ref_len) {
            // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to be reference length
            fragment_len = ref_len;
        } else {
            fragment_len = 0;
            while (fragment_len < _art_params.read_len || fragment_len > _ref_seq.length()) {
                fragment_len = static_cast<int>(gaussian(rng));
            }
        }
    }
    auto pos_1 = is_amplicon ? 0 : static_cast<int>(floor((ref_len - fragment_len) * r_prob()));
    auto pos_2 = is_amplicon ? 0 : ref_len - pos_1 - fragment_len;
    int slen_1 = read_1.generate_indels(_art_params.read_len, true);
    int slen_2 = read_2.generate_indels(_art_params.read_len, false);

    // ensure get a fixed read length
    if ((pos_1 + _art_params.read_len - slen_1) > ref_len) {
        slen_1 = read_1.generate_indels_2(_art_params.read_len, true);
    }
    if ((pos_2 + _art_params.read_len - slen_2) > ref_len) {
        slen_2 = read_2.generate_indels_2(_art_params.read_len, false);
    }

    bool is_plus_strand = is_amplicon || r_prob() < 0.5;

    // string read_info_1;
    // string read_info_2;
    // int aln_start;
    // int aln_end;
    if (is_plus_strand) {
        //    |----------->
        // ------------------------------------
        //                   <-----------|
        read_1.is_plus_strand = true;
        read_1.seq_ref = _ref_seq.substr(pos_1, _art_params.read_len - slen_1);
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % pos_1 % (_art_params.read_len - slen_1 + pos_1) % "+").str();
        read_2.is_plus_strand = false;
        read_2.seq_ref = _ref_seq_cmp.substr(pos_2, _art_params.read_len - slen_2);
        // read_info_2 = (boost::format("%s:%d-%d:%s") % _id % (ref_len - pos_2 + 1) % (ref_len - pos_2 + 1 - _art_params.read_len + slen_2) % "-").str();
        // aln_start = pos_1;
        // aln_end = ref_len - pos_2 + 1;
    } else {
        //                   <-----------|
        // ------------------------------------
        //    |----------->
        read_1.is_plus_strand = false;
        read_1.seq_ref = _ref_seq_cmp.substr(pos_1, _art_params.read_len - slen_1);
        // read_info_1 = (boost::format("%s:%d-%d:%s") % _id % (ref_len - pos_1 + 1) % (ref_len - pos_1 + 1 - _art_params.read_len + slen_1) % "-").str();
        read_2.is_plus_strand = true;
        read_2.seq_ref = _ref_seq.substr(pos_2, _art_params.read_len - slen_2);
        // read_info_2 = (boost::format("%s:%d-%d:%s") % _id % pos_2 % (_art_params.read_len - slen_2 + pos_2) % "+").str();
        // aln_start = pos_2;
        // aln_end = ref_len - pos_1 + 1;
    }
    // std::cout << "R1: " << read_info_1 << ", R2: " << read_info_2 << " Dist: " << aln_end - aln_start + 1 << endl;
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
    _ref_seq_cmp.resize(ref_len);
    boost::algorithm::to_upper(ref_seq);
    for (auto i = 0; i < ref_len; i++) {
        auto k = ref_len - i - 1; // Reverse complementary position
        switch (ref_seq[i]) {
        case 'A':
            _ref_seq_cmp[k] = 'T';
            break;
        case 'C':
            _ref_seq_cmp[k] = 'G';
            break;
        case 'G':
            _ref_seq_cmp[k] = 'C';
            break;
        case 'T':
            _ref_seq_cmp[k] = 'A';
            break;
        default:
            _ref_seq_cmp[k] = 'N';
        }
    }
    _valid_region = ref_len - _art_params.read_len;
}
