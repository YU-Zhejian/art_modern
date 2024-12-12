#include "Empdist.hh"
#include "ArtConstants.hh"
#include "random_generator.hh"
#include "utils/mpi_utils.hh"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <fstream>
#include <sstream>

namespace labw::art_modern {

void shift_emp(Empdist::dist_type& map_to_process, const int q_shift, const int min_qual, const int max_qual)
{
    for (auto& i : map_to_process) {
        for (auto& j : i) {
            j.second += q_shift;
            j.second = std::min(std::max(j.second, min_qual), max_qual);
        }
    }
}

Empdist::Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual)
    : sep_qual_(sep_qual)
{
    read_emp_dist_(emp_filename_1, true);
    if (!emp_filename_2.empty()) {
        read_emp_dist_(emp_filename_2, false);
    }
    if (sep_qual_) {
        if (a_qual_dist_first.size() != g_qual_dist_first.size() || g_qual_dist_first.size() != c_qual_dist_first.size()
            || c_qual_dist_first.size() != t_qual_dist_first.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
            abort_mpi();
        }
        if (a_qual_dist_second.size() != g_qual_dist_second.size()
            || g_qual_dist_second.size() != c_qual_dist_second.size()
            || c_qual_dist_second.size() != t_qual_dist_second.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
            abort_mpi();
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Read quality profile loaded successfully.";
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

// generate quality vector from dist of one read from pair-end [default first
// read]
void Empdist::get_read_qual(std::vector<int>& qual, int len, Rprob& rprob, const bool first) const
{
    qual.resize(len);
    const auto& qual_dist = first ? qual_dist_first : qual_dist_second;
    rprob.rand_quality(qual);
    for (auto i = 0; i < len; i++) {
        qual[i] = qual_dist[i].lower_bound(qual[i])->second;
    }
}

void Empdist::get_read_qual_sep_1(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const
{
    const auto len = seq.size();
    qual.resize(len);

    if (a_qual_dist_first.size() < len || t_qual_dist_first.size() < len || g_qual_dist_first.size() < len
        || c_qual_dist_first.size() < len) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << len
                                 << ") exceeds the length of the read quality profile (A: " << a_qual_dist_first.size()
                                 << ", C: " << c_qual_dist_first.size() << ", G: " << g_qual_dist_first.size()
                                 << ", T: " << t_qual_dist_first.size()

                                 << ")";
        abort_mpi();
    }

    rprob.rand_quality(qual);

    for (decltype(seq.size()) i = 0; i < len; i++) {
        if (seq[i] == 'A') {
            qual[i] = a_qual_dist_first[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'C') {
            qual[i] = c_qual_dist_first[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'G') {
            qual[i] = g_qual_dist_first[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'T') {
            qual[i] = t_qual_dist_first[i].lower_bound(qual[i])->second;
        } else {
            // return random quality less than 10
            qual[i] = rprob.rand_quality_less_than_10();
        }
    }
}

void Empdist::get_read_qual_sep_2(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const
{
    const auto len = seq.size();
    qual.resize(len);

    if (a_qual_dist_second.size() < len || t_qual_dist_second.size() < len || g_qual_dist_second.size() < len
        || c_qual_dist_second.size() < len) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 2nd read (" << len
                                 << ") exceeds the "
                                    "length of the read quality profile (A: "
                                 << a_qual_dist_second.size() << ", C: " << c_qual_dist_second.size()
                                 << ", G: " << g_qual_dist_second.size() << ", T: " << t_qual_dist_second.size() << ")";
        abort_mpi();
    }
    rprob.rand_quality(qual);
    for (size_t i = 0; i < len; i++) {
        if (seq[i] == 'A') {
            qual[i] = a_qual_dist_second[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'C') {
            qual[i] = c_qual_dist_second[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'G') {
            qual[i] = g_qual_dist_second[i].lower_bound(qual[i])->second;
        } else if (seq[i] == 'T') {
            qual[i] = t_qual_dist_second[i].lower_bound(qual[i])->second;
        } else {
            // return random quality less than 10
            qual[i] = rprob.rand_quality_less_than_10();
        }
    }
}

void Empdist::read_emp_dist_(std::istream& input, const bool is_first)
{
    int linenum = 0;
    int read_pos;
    char alt_read_pos;
    char leading_base;
    long t_uint;
    std::string line;
    int t_int;
    std::vector<int> qual;
    std::map<int, int, std::less<>> dist;
    std::vector<long> count;

    while (!input.eof()) {
        std::getline(input, line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        leading_base = line[0];
        if (!((sep_qual_ && ART_ACGT_STR.find(leading_base) != std::string::npos)
                || (!sep_qual_ && leading_base == '.'))) {
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

        auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;
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

}