#include "art/lib/Empdist.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/BuiltinProfile.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/utils/mpi_utils.hh"

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
    void shift_emp(Empdist::dist_type& map_to_process, const int q_shift, const int min_qual, const int max_qual)
    {
        for (auto& i : map_to_process) {
            for (auto& [fst, snd] : i) {
                snd += q_shift;
                snd = std::min(std::max(snd, min_qual), max_qual);
            }
        }
    }
} // namespace

Empdist::Empdist(const BuiltinProfile& builtin_profile, bool sep_qual, bool is_pe, std::size_t read_len)
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
    print_();
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
    print_();
}

// generate quality vector from dist of one read from pair-end [default first
// read]
void Empdist::get_read_qual(std::vector<am_qual_t>& qual, const int len, Rprob& rprob, const bool first) const
{
    const auto& qual_dist = first ? qual_dist_first : qual_dist_second;
    rprob.rand_quality_dist();
    for (auto i = 0; i < len; i++) {
        // TODO: This line of code have catastrophic locality.
        qual[i] = qual_dist[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
    }
}

void Empdist::get_read_qual_sep_1(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const
{
    const auto len = seq.size();

    rprob.rand_quality_dist();

    for (decltype(seq.size()) i = 0; i < len; i++) {
        switch (seq[i]) {
        case 'A':
            qual[i] = a_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'C':
            qual[i] = c_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'G':
            qual[i] = g_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'T':
            qual[i] = t_qual_dist_first[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
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
            qual[i] = a_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'C':
            qual[i] = c_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'G':
            qual[i] = g_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        case 'T':
            qual[i] = t_qual_dist_second[i].lower_bound(rprob.tmp_qual_dists_[i])->second;
            break;
        default:
            qual[i] = rprob.rand_quality_less_than_10();
        }
    }
}

void Empdist::read_emp_dist_(std::istream& input, const bool is_first)
{
    int linenum = 0;
    int read_pos = 0;
    char alt_read_pos = 0;
    char leading_base = 0;
    am_qual_dist_t t_uint = 0;
    std::string line;
    int t_int = 0;
    std::vector<int> qual;
    dist_map_type dist;
    std::vector<am_qual_dist_t> count;
    int qmin = std::numeric_limits<int>::max();
    int qmax = std::numeric_limits<int>::min();

    while (!input.eof()) {
        std::getline(input, line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        leading_base = line[0];
        if ((!sep_qual_ || ART_ACGT_STR.find(leading_base) == std::string::npos)
            && (sep_qual_ || leading_base != '.')) {
            continue;
        }

        std::istringstream ss(line);
        ss >> alt_read_pos;
        ss >> read_pos;

        if (read_pos != linenum) {
            linenum = 0;
            if (read_pos != linenum) {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal error (1): Wrong format of input distribution.";
                abort_mpi();
            }
        }
        qual.clear();
        while (ss >> t_int) {
            qual.emplace_back(t_int);
        }
        qmin = std::min(qmin, *std::min_element(qual.begin(), qual.end()));
        qmax = std::max(qmax, *std::min_element(qual.begin(), qual.end()));

        std::getline(input, line);
        ss.clear();
        ss.str(line);
        ss >> alt_read_pos;
        ss >> read_pos;

        if (read_pos != linenum) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal error (2): Wrong format of input distribution.";
            abort_mpi();
        }

        count.clear();
        count.reserve(qual.size());

        while (ss >> t_uint) {
            count.emplace_back(t_uint);
        }

        if (count.size() != qual.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal error (3): Wrong format of input distribution. count.size = "
                                     << count.size() << "; qual_.size = " << qual.size();
            abort_mpi();
        }

        const auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;
        dist.clear();

        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
        }
        if (!dist.empty()) {
            linenum++;
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
    }

    if (linenum == 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Profile empty!";
        abort_mpi();
    }
    BOOST_LOG_TRIVIAL(info) << "QRange for R" << (is_first ? 1 : 2) << ": [" << qmin << ", " << qmax << "].";
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

void Empdist::shift_all_emp(
    const bool sep_flag, const int q_shift_1, const int q_shift_2, const int min_qual, const int max_qual)
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
            BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: The length of 1st read in each qual dist is not equal!";
            print_();
            abort_mpi();
        }
        if (a_qual_dist_second.size() != g_qual_dist_second.size()
            || g_qual_dist_second.size() != c_qual_dist_second.size()
            || c_qual_dist_second.size() != t_qual_dist_second.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: The length of 2nd read in each qual dist is not equal!";
            print_();
            abort_mpi();
        }
    }
    if (sep_qual_) {
        if (a_qual_dist_first.size() < read_len_) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << read_len_
                                     << ") exceeds the "
                                        "length of the read quality profile.";
            print_();
            abort_mpi();
        }
        if (is_pe_) {
            if (a_qual_dist_second.size() < read_len_) {
                BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 2nd read (" << read_len_
                                         << ") exceeds the "
                                            "length of the read quality profile.";
                print_();
                abort_mpi();
            }
        }
    } else {
        if (qual_dist_first.size() < read_len_) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << read_len_
                                     << ") exceeds the "
                                        "length of the read quality profile.";
            print_();
            abort_mpi();
        }
        if (is_pe_) {
            if (qual_dist_second.size() < read_len_) {
                BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 2nd read (" << read_len_
                                         << ") exceeds the "
                                            "length of the read quality profile ("
                                         << qual_dist_second.size() << ")";
                print_();
                abort_mpi();
            }
        }
    }
}

void Empdist::print_() const
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

} // namespace labw::art_modern
