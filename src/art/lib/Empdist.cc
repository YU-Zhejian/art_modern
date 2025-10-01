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

#include "art/lib/Empdist.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/BuiltinProfile.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/arithmetic_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace labw::art_modern {

namespace {
    void shift_emp(
        Empdist::dist_type& map_to_process, const am_qual_t q_shift, const am_qual_t min_qual, const am_qual_t max_qual)
    {
        for (auto& i : map_to_process) {
            for (auto& [fst, snd] : i) {
                snd += q_shift;
                snd = am_min(am_max(snd, min_qual), max_qual);
            }
        }
    }
} // namespace

Empdist::Empdist(
    const BuiltinProfile& builtin_profile, const bool sep_qual, const bool is_pe, const std::size_t read_len)
    : sep_qual_(sep_qual)
    , is_pe_(is_pe)
    , read_len_(read_len)
{
    std::istringstream ss(builtin_profile.r1_profile);
    read_emp_dist_(ss, true);
    if (!builtin_profile.r2_profile.empty()) {
        std::istringstream ss2(builtin_profile.r2_profile);
        read_emp_dist_(ss2, false);
    }
    validate_();
    BOOST_LOG_TRIVIAL(info) << "Read quality profile loaded successfully.";
    log();
}

Empdist::Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, const bool sep_qual,
    const bool is_pe, const std::size_t read_len)
    : sep_qual_(sep_qual)
    , is_pe_(is_pe)
    , read_len_(read_len)
{
    read_emp_dist_(emp_filename_1, true);
    if (!emp_filename_2.empty()) {
        read_emp_dist_(emp_filename_2, false);
    }
    validate_();
    BOOST_LOG_TRIVIAL(info) << "Read quality profile loaded successfully.";
    log();
}

// generate quality vector from dist of one read from pair-end [default first
// read]
void Empdist::get_read_qual(std::vector<am_qual_t>& qual, Rprob& rprob, const bool first) const
{
#ifdef USE_WALKER_QUALGEN
    const auto& qual_dist_idx = first ? qual_dist_first_idx : qual_dist_second_idx;
    rprob.r_probs();
    for (std::size_t i = 0; i < read_len_; i++) {
        // TODO: This line of code have catastrophic locality.
        qual[i] = qual_dist_idx[i].gen_qual(rprob.tmp_probs_[i]);
    }
#else
    const auto& qual_dist = first ? qual_dist_first : qual_dist_second;
    rprob.rand_quality_dist();
    for (std::size_t i = 0; i < read_len_; i++) {
        // TODO: This line of code have catastrophic locality.
        qual[i] = qual_dist[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
    }
#endif
}

void Empdist::get_read_qual_sep_1(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const
{
    const auto len = seq.size();

#ifdef USE_WALKER_QUALGEN
    rprob.r_probs();
#else
    rprob.rand_quality_dist();
#endif
    for (decltype(seq.size()) i = 0; i < len; i++) {
        switch (seq[i]) {
        case 'A':
#ifdef USE_WALKER_QUALGEN
            qual[i] = a_qual_dist_first_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = a_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'C':
#ifdef USE_WALKER_QUALGEN
            qual[i] = c_qual_dist_first_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = c_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'G':
#ifdef USE_WALKER_QUALGEN
            qual[i] = g_qual_dist_first_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = g_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'T':
#ifdef USE_WALKER_QUALGEN
            qual[i] = t_qual_dist_first_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = t_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        default:
            qual[i] = rprob.rand_quality_less_than_10();
        }
    }
}

void Empdist::get_read_qual_sep_2(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const
{
    const auto len = seq.size();

    rprob.rand_quality_dist();
    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
        case 'A':
#ifdef USE_WALKER_QUALGEN
            qual[i] = a_qual_dist_second_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = a_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'C':
#ifdef USE_WALKER_QUALGEN
            qual[i] = c_qual_dist_second_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = c_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'G':
#ifdef USE_WALKER_QUALGEN
            qual[i] = g_qual_dist_second_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = g_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        case 'T':
#ifdef USE_WALKER_QUALGEN
            qual[i] = t_qual_dist_second_idx[i].gen_qual(rprob.tmp_probs_[i]);
#else
            qual[i] = t_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
#endif
            break;
        default:
            qual[i] = rprob.rand_quality_less_than_10();
        }
    }
}

void Empdist::read_emp_dist_(std::istream& input, const bool is_first)
{
    std::size_t actual_line_no = 0;
    std::size_t n_lines_parsed = 0;
    /** where we are on the read **/
    std::size_t read_pos = 0;
    char leading_base = 0;
    std::string line;

    // ChatGPT: When using stringstream >> int8_t,
    // the extraction operator treats int8_t as a char,
    // so it reads a single character, not a number.
    // This can lead to unexpected results.
    int tmp_qual = 0;
    am_qual_count_t tmp_count = 0;

    std::vector<am_qual_t> qual;
    std::vector<am_qual_count_t> count;

    dist_map_type dist;
    am_qual_t qmin = std::numeric_limits<am_qual_t>::max();
    am_qual_t qmax = std::numeric_limits<am_qual_t>::min();

    while (!input.eof()) {
        std::getline(input, line);
        actual_line_no++;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        leading_base = line[0];
        if ((sep_qual_ && (leading_base == 'A' || leading_base == 'C' || leading_base == 'G' || leading_base == 'T'))
            || (!sep_qual_ && leading_base == '.')) {
            // Normal condition.
            // Note that the current version does not support N now
        } else {
            continue;
        }

        std::istringstream ss(line);

        ss >> leading_base;
        ss >> read_pos;

        if (read_pos != n_lines_parsed) {
            n_lines_parsed = 0;
            if (read_pos != n_lines_parsed) {
                BOOST_LOG_TRIVIAL(fatal) << "R" << (is_first ? 1 : 2) << "L" << std::to_string(actual_line_no)
                                         << ": Fatal error (1): Wrong format of input distribution.";
                BOOST_LOG_TRIVIAL(fatal) << "line=" << line;
                abort_mpi();
            }
        }

        qual.clear();
        while (ss >> tmp_qual) {
            qual.emplace_back(tmp_qual);
        }

        qmin = am_min(qmin, *std::min_element(qual.begin(), qual.end()));
        qmax = am_max(qmax, *std::max_element(qual.begin(), qual.end()));

        std::getline(input, line);
        actual_line_no++;

        ss.clear();
        ss.str(line);
        ss >> leading_base;
        ss >> read_pos;

        if (read_pos != n_lines_parsed) {
            BOOST_LOG_TRIVIAL(fatal) << "R" << (is_first ? 1 : 2) << "L" << std::to_string(actual_line_no)
                                     << ": Fatal error (2): Wrong format of input distribution.";
            BOOST_LOG_TRIVIAL(fatal) << "read_pos=" << read_pos << "; n_lines_parsed=" << n_lines_parsed;
            BOOST_LOG_TRIVIAL(fatal) << "line=" << line;
            abort_mpi();
        }

        count.clear();
        count.reserve(qual.size());

        while (ss >> tmp_count) {
            count.emplace_back(tmp_count);
        }

        if (count.size() != qual.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "R" << (is_first ? 1 : 2) << "L" << std::to_string(actual_line_no)
                                     << ": Fatal error (3): Wrong format of input distribution. count.size = "
                                     << count.size() << "; qual.size = " << qual.size();
            BOOST_LOG_TRIVIAL(fatal) << "qual=" << vec2str(qual);
            BOOST_LOG_TRIVIAL(fatal) << "count=" << vec2str(count);
            abort_mpi();
        }

        const auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;
        dist.clear();

        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
        }
        n_lines_parsed++;
        if (!sep_qual_) {
            (is_first ? qual_dist_first : qual_dist_second).emplace_back(dist);
        } else if (leading_base == 'A') {
            (is_first ? a_qual_dist_first : a_qual_dist_second).emplace_back(dist);
        } else if (leading_base == 'T') {
            (is_first ? t_qual_dist_first : t_qual_dist_second).emplace_back(dist);
        } else if (leading_base == 'G') {
            (is_first ? g_qual_dist_first : g_qual_dist_second).emplace_back(dist);
        } else if (leading_base == 'C') {
            (is_first ? c_qual_dist_first : c_qual_dist_second).emplace_back(dist);
        } else {
            BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
            abort_mpi();
        }
    }

    if (n_lines_parsed == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Profile empty!";
        abort_mpi();
    }
    BOOST_LOG_TRIVIAL(info) << "QRange for R" << (is_first ? 1 : 2) << ": [" << std::to_string(qmin) << ", "
                            << std::to_string(qmax) << "].";
}

void Empdist::read_emp_dist_(const std::string& infile, const bool is_first)
{
    std::ifstream distss(infile);
    if (!distss) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Cannot open the distribution file: '" << infile << "'";
        abort_mpi();
    }
    read_emp_dist_(distss, is_first);
}

void Empdist::shift_all_emp(const bool sep_flag, const am_qual_t q_shift_1, const am_qual_t q_shift_2,
    const am_qual_t min_qual, const am_qual_t max_qual)
{
    if (!sep_flag) {
        shift_emp(qual_dist_first, q_shift_1, min_qual, max_qual);
        shift_emp(qual_dist_second, q_shift_2, min_qual, max_qual);
    } else {
        shift_emp(a_qual_dist_first, q_shift_1, min_qual, max_qual);
        shift_emp(a_qual_dist_second, q_shift_2, min_qual, max_qual);
        shift_emp(t_qual_dist_first, q_shift_1, min_qual, max_qual);
        shift_emp(t_qual_dist_second, q_shift_2, min_qual, max_qual);
        shift_emp(c_qual_dist_first, q_shift_1, min_qual, max_qual);
        shift_emp(c_qual_dist_second, q_shift_2, min_qual, max_qual);
        shift_emp(g_qual_dist_first, q_shift_1, min_qual, max_qual);
        shift_emp(g_qual_dist_second, q_shift_2, min_qual, max_qual);
    }
}

void Empdist::validate_() const
{
    if (sep_qual_) {
        if (a_qual_dist_first.size() != g_qual_dist_first.size() || g_qual_dist_first.size() != c_qual_dist_first.size()
            || c_qual_dist_first.size() != t_qual_dist_first.size()) {
            BOOST_LOG_TRIVIAL(warning) << "The length of 1st read in each qual dist is not equal!";
        }
        if (a_qual_dist_second.size() != g_qual_dist_second.size()
            || g_qual_dist_second.size() != c_qual_dist_second.size()
            || c_qual_dist_second.size() != t_qual_dist_second.size()) {
            BOOST_LOG_TRIVIAL(warning) << "The length of 2nd read in each qual dist is not equal!";
        }
    }
    if (sep_qual_) {
        if (a_qual_dist_first.size() < read_len_) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << read_len_
                                     << ") exceeds the "
                                        "length of the read quality profile.";
            log();
            abort_mpi();
        }
        if (is_pe_) {
            if (a_qual_dist_second.size() < read_len_) {
                BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 2nd read (" << read_len_
                                         << ") exceeds the "
                                            "length of the read quality profile.";
                log();
                abort_mpi();
            }
        }
    } else {
        if (qual_dist_first.size() < read_len_) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << read_len_
                                     << ") exceeds the "
                                        "length of the read quality profile.";
            log();
            abort_mpi();
        }
        if (is_pe_) {
            if (qual_dist_second.size() < read_len_) {
                BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 2nd read (" << read_len_
                                         << ") exceeds the "
                                            "length of the read quality profile ("
                                         << qual_dist_second.size() << ")";
                log();
                abort_mpi();
            }
        }
    }
}

void Empdist::log() const
{
    if (sep_qual_) {
        BOOST_LOG_TRIVIAL(info) << "Read quality profile size for R1: A: " << a_qual_dist_first.size()
                                << ", C: " << c_qual_dist_first.size() << ", G: " << g_qual_dist_first.size()
                                << ", T: " << t_qual_dist_first.size();
        BOOST_LOG_TRIVIAL(info) << "Read quality profile size for R2: A: " << a_qual_dist_second.size()
                                << ", C: " << c_qual_dist_second.size() << ", G: " << g_qual_dist_second.size()
                                << ", T: " << t_qual_dist_second.size();
    } else {
        BOOST_LOG_TRIVIAL(info) << "Read quality profile size for R1: " << qual_dist_first.size();
        BOOST_LOG_TRIVIAL(info) << "Read quality profile size for R2: " << qual_dist_second.size();
    }
}

void Empdist::index()
{
#ifdef USE_WALKER_QUALGEN
    if (sep_qual_) {
        for (const auto& dist : a_qual_dist_first) {
            a_qual_dist_first_idx.emplace_back(dist);
        }
        for (const auto& dist : a_qual_dist_second) {
            a_qual_dist_second_idx.emplace_back(dist);
        }
        for (const auto& dist : c_qual_dist_first) {
            c_qual_dist_first_idx.emplace_back(dist);
        }
        for (const auto& dist : c_qual_dist_second) {
            c_qual_dist_second_idx.emplace_back(dist);
        }
        for (const auto& dist : g_qual_dist_first) {
            g_qual_dist_first_idx.emplace_back(dist);
        }
        for (const auto& dist : g_qual_dist_second) {
            g_qual_dist_second_idx.emplace_back(dist);
        }
        for (const auto& dist : t_qual_dist_first) {
            t_qual_dist_first_idx.emplace_back(dist);
        }
        for (const auto& dist : t_qual_dist_second) {
            t_qual_dist_second_idx.emplace_back(dist);
        }
    } else {

        for (const auto& dist : qual_dist_first) {
            qual_dist_first_idx.emplace_back(dist);
        }
        for (const auto& dist : qual_dist_second) {
            qual_dist_second_idx.emplace_back(dist);
        }
    }

#endif
}

} // namespace labw::art_modern
