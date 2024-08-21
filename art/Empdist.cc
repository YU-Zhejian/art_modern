#include "Empdist.hh"
#include "ArtConstants.hh"
#include "random_generator.hh"
#include <boost/log/trivial.hpp>
#include <fstream>
#include <sstream>

using namespace std;
namespace labw {
namespace art_modern {
    Empdist::Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual)
        : sep_qual_(sep_qual)
    {
        read_emp_dist_(emp_filename_1, true);
        if (!emp_filename_2.empty()) {
            read_emp_dist_(emp_filename_2, false);
        }
    }

    // generate quality vector from dist of one read from pair-end [default first
    // read]
    void Empdist::get_read_qual(std::vector<int>& qual, int len, Rprob& rprob, bool first) const
    {
        qual.resize(len);
        auto qual_dist = first ? qual_dist_first : qual_dist_second;
        int cumCC;
        map<int, int>::iterator it;
        for (int i = 0; i < len; i++) {
            cumCC = rprob.rand_quality();
            it = qual_dist[i].lower_bound(cumCC);
            qual[i] = it->second;
        }
    }

    void Empdist::get_read_qual_sep_1(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const
    {
        const auto len = seq.size();
        qual.resize(len);

        if (a_qual_dist_first.size() < len || t_qual_dist_first.size() < len || g_qual_dist_first.size() < len
            || c_qual_dist_first.size() < len) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length exceeds the "
                                        "length of the read quality profile";
            exit(EXIT_FAILURE);
        }

        int cumCC;

        for (size_t i = 0; i < len; i++) {
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
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length exceeds the "
                                        "length of the read quality profile";
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

    void Empdist::read_emp_dist_(istream& input, bool is_first)
    {
        int linenum = 0;
        int read_pos;
        char alt_read_pos;

        while (!input.eof()) {
            bool a_flag = false;
            bool t_flag = false;
            bool g_flag = false;
            bool c_flag = false;

            string aLine;
            getline(input, aLine);
            if (aLine.empty())
                continue;
            if (aLine[0] == '#') {
                continue;
            }

            if (!sep_qual_ && (aLine[0] != '.')) {
                continue;
            }
            if (sep_qual_) {
                if (aLine[0] == 'A') {
                    a_flag = true;
                } else if (aLine[0] == 'T') {
                    t_flag = true;
                } else if (aLine[0] == 'G') {
                    g_flag = true;
                } else if (aLine[0] == 'C') {
                    c_flag = true;
                } else if (aLine[0] == 'N' || aLine[0] == '.') {
                    continue;
                }
            }

            istringstream ss(aLine);
            ss >> alt_read_pos;
            ss >> read_pos;

            if (read_pos != linenum) {
                linenum = 0;
                if (read_pos != linenum) {
                    BOOST_LOG_TRIVIAL(fatal) << "Fatal error (1): Wrong format of input distribution.";
                    exit(EXIT_FAILURE);
                }
            }

            int t_int;
            vector<int> qual;

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

            long t_uint;
            vector<long> count;
            count.reserve(qual.size());

            while (ss >> t_uint) {
                count.emplace_back(t_uint);
            }

            if (count.size() != qual.size()) {
                BOOST_LOG_TRIVIAL(fatal) << "Fatal error (3): Wrong format of input distribution. count.size = "
                                         << count.size() << "; qual.size = " << qual.size();
                exit(EXIT_FAILURE);
            }

            double denom = static_cast<double>(count[count.size() - 1]) / MAX_DIST_NUMBER;
            map<int, int> dist;

            for (size_t i = 0; i < count.size(); i++) {
                auto cc = static_cast<int>(ceil(static_cast<double>(count[i]) / denom));
                dist[cc] = qual[i];
            }
            if (!dist.empty()) {
                linenum++;
                if (!sep_qual_ && is_first) {
                    qual_dist_first.emplace_back(dist);
                } else if (!sep_qual_) {
                    qual_dist_second.emplace_back(dist);
                } else if (a_flag && is_first) {
                    a_qual_dist_first.emplace_back(dist);
                } else if (a_flag) {
                    a_qual_dist_second.emplace_back(dist);
                } else if (t_flag && is_first) {
                    t_qual_dist_first.emplace_back(dist);
                } else if (t_flag) {
                    t_qual_dist_second.emplace_back(dist);
                } else if (g_flag && is_first) {
                    g_qual_dist_first.emplace_back(dist);
                } else if (g_flag) {
                    g_qual_dist_second.emplace_back(dist);
                } else if (c_flag && is_first) {
                    c_qual_dist_first.emplace_back(dist);
                } else if (c_flag) {
                    c_qual_dist_second.emplace_back(dist);
                } else {
                    BOOST_LOG_TRIVIAL(fatal) << "Unexpected Error: Profile was not read in correctly.";
                    exit(EXIT_FAILURE);
                }
            }
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

        if (linenum == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Profile empty!";
            exit(EXIT_FAILURE);
        }
    }

    void Empdist::read_emp_dist_(const std::string& infile, bool is_first)
    {
        ifstream distss(infile.c_str());
        if (!distss) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Cannot open the distribution file: '" << infile << "'";
            exit(EXIT_FAILURE);
        }
        read_emp_dist_(distss, is_first);
    }

}
}