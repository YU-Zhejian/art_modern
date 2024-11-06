#include "Empdist.hh"
#include "ArtConstants.hh"
#include "random_generator.hh"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <fstream>
#include <sstream>

namespace labw {
namespace art_modern {
    Empdist::Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual)
        : sep_qual_(sep_qual)
    {
        read_emp_dist_(emp_filename_1, true);
        if (!emp_filename_2.empty()) {
            read_emp_dist_(emp_filename_2, false);
        }
        if (sep_qual_) {
            if (a_qual_dist_first.size() != g_qual_dist_first.size()
                || g_qual_dist_first.size() != c_qual_dist_first.size()
                || c_qual_dist_first.size() != t_qual_dist_first.size()) {
                BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
                exit(EXIT_FAILURE);
            }
            if (a_qual_dist_second.size() != g_qual_dist_second.size()
                || g_qual_dist_second.size() != c_qual_dist_second.size()
                || c_qual_dist_second.size() != t_qual_dist_second.size()) {
                BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
                exit(EXIT_FAILURE);
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
    void Empdist::get_read_qual(std::vector<int>& qual, int len, Rprob& rprob, bool first) const
    {
        qual.resize(len);
        const auto& qual_dist = first ? qual_dist_first : qual_dist_second;
        int cumCC;
        for (int i = 0; i < len; i++) {
            cumCC = rprob.rand_quality();
            qual[i] = qual_dist[i].lower_bound(cumCC)->second;
        }
    }

    void Empdist::get_read_qual_sep_1(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const
    {
        const auto len = seq.size();
        qual.resize(len);

        if (a_qual_dist_first.size() < len || t_qual_dist_first.size() < len || g_qual_dist_first.size() < len
            || c_qual_dist_first.size() < len) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length of 1st read (" << len
                                     << ") exceeds the length of the read quality profile (A: "
                                     << a_qual_dist_first.size() << ", C: " << c_qual_dist_first.size()
                                     << ", G: " << g_qual_dist_first.size() << ", T: " << t_qual_dist_first.size()

                                     << ")";
            exit(EXIT_FAILURE);
        }

        int cumCC;

        for (auto i = 0; i < len; i++) {
            cumCC = rprob.rand_quality();
            if (seq[i] == 'A') {
                auto it = a_qual_dist_first[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'C') {
                auto it = c_qual_dist_first[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'G') {
                auto it = g_qual_dist_first[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'T') {
                auto it = t_qual_dist_first[i].lower_bound(cumCC);
                qual[i] = it->second;
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
                                     << ", G: " << g_qual_dist_second.size() << ", T: " << t_qual_dist_second.size()
                                     << ")";
            exit(EXIT_FAILURE);
        }
        int cumCC;
        for (size_t i = 0; i < len; i++) {
            cumCC = rprob.rand_quality(); // (int)ceil(rprob.r_prob() * max_dist_number);
            if (seq[i] == 'A') {
                auto it = a_qual_dist_second[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'C') {
                auto it = c_qual_dist_second[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'G') {
                auto it = g_qual_dist_second[i].lower_bound(cumCC);
                qual[i] = it->second;
            } else if (seq[i] == 'T') {
                auto it = t_qual_dist_second[i].lower_bound(cumCC);
                qual[i] = it->second;
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
        bool a_flag;
        bool t_flag;
        bool g_flag;
        bool c_flag;
        long t_uint;
        std::string aLine;

        int t_int;
        std::vector<int> qual;
        std::map<int, int> dist;
        std::vector<long> count;

        while (!input.eof()) {
            a_flag = false;
            t_flag = false;
            g_flag = false;
            c_flag = false;

            std::getline(input, aLine);
            if (aLine.empty() || aLine[0] == '#') {
                continue;
            }
            if (sep_qual_) {
                switch (aLine[0]) {
                case 'A':
                    a_flag = true;
                    break;
                case 'T':
                    t_flag = true;
                    break;
                case 'G':
                    g_flag = true;
                    break;
                case 'C':
                    c_flag = true;
                    break;
                default:
                    continue;
                }
            } else {
                if (aLine[0] != '.') {
                    continue;
                }
            }

            std::istringstream ss(aLine);
            ss >> alt_read_pos;
            ss >> read_pos;

            if (read_pos != linenum) {
                linenum = 0;
                if (read_pos != linenum) {
                    BOOST_LOG_TRIVIAL(fatal) << "Fatal error (1): Wrong format of input distribution.";
                    exit(EXIT_FAILURE);
                }
            }
            qual.clear();
            while (ss >> t_int) {
                qual.emplace_back(t_int);
            }

            getline(input, aLine);
            ss.clear();
            ss.str(aLine);
            ss >> alt_read_pos;
            ss >> read_pos;

            if (read_pos != linenum) {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal error (2): Wrong format of input distribution.";
                exit(EXIT_FAILURE);
            }

            count.clear();
            count.reserve(qual.size());

            while (ss >> t_uint) {
                count.emplace_back(t_uint);
            }

            if (count.size() != qual.size()) {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal error (3): Wrong format of input distribution. count.size = "
                                         << count.size() << "; qual.size = " << qual.size();
                exit(EXIT_FAILURE);
            }

            auto denom = static_cast<double>(count[count.size() - 1]) / MAX_DIST_NUMBER;
            dist.clear();

            for (auto i = 0; i < count.size(); i++) {
                dist[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
            }
            if (!dist.empty()) {
                linenum++;
                if (!sep_qual_) {
                    (is_first ? qual_dist_first : qual_dist_second).emplace_back(dist);
                } else if (a_flag) {
                    (is_first ? a_qual_dist_first : a_qual_dist_second).emplace_back(dist);
                } else if (t_flag) {
                    (is_first ? t_qual_dist_first : t_qual_dist_second).emplace_back(dist);
                } else if (g_flag) {
                    (is_first ? g_qual_dist_first : g_qual_dist_second).emplace_back(dist);
                } else if (c_flag) {
                    (is_first ? c_qual_dist_first : c_qual_dist_second).emplace_back(dist);
                } else {
                    BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
                    exit(EXIT_FAILURE);
                }
            }
        }

        if (linenum == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Profile empty!";
            exit(EXIT_FAILURE);
        }
    }

    void Empdist::read_emp_dist_(const std::string& infile, const bool is_first)
    {
        std::ifstream distss(infile);
        if (!distss) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Cannot open the distribution file: '" << infile << "'";
            exit(EXIT_FAILURE);
        }
        read_emp_dist_(distss, is_first);
    }

}
}