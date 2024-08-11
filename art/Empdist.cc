#include "Empdist.hh"
#include <boost/log/trivial.hpp>
#include <fstream>
#include <sstream>

using namespace std;
namespace labw {
namespace art_modern {

    Empdist::Empdist(const std::string& emp_filename_1,
        const std::string& emp_filename_2, bool sep_qual)
        : _sep_qual(sep_qual)
    {
        read_emp_dist(emp_filename_1, true);
        if (!emp_filename_2.empty()) {
            read_emp_dist(emp_filename_2, false);
        }
    }

    // generate quali vector from dist of one read from pair-end [default first
    // read]
    std::vector<int> Empdist::get_read_qual(int len, bool first) const
    {
        std::vector<int> read_qual;
        auto qual_dist = first ? qual_dist_first : qual_dist_second;
        int cumCC;
        map<int, int>::iterator it;
        for (int i = 0; i < len; i++) {
            cumCC = (int)ceil(r_prob() * max_dist_number);
            it = qual_dist[i].lower_bound(cumCC);
            read_qual.push_back(it->second);
        }
        return read_qual;
    }

    vector<int> Empdist::get_read_qual_sep_1(const string& seq) const
    {
        vector<int> read_qual;
        const auto len = static_cast<int>(seq.size());
        if (a_qual_dist_first.size() < len || t_qual_dist_first.size() < len || g_qual_dist_first.size() < len
            || c_qual_dist_first.size() < len) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length exceeds the "
                                        "length of the read quality profile";
            exit(EXIT_FAILURE);
        }

        int cumCC;

        for (int i = 0; i < len; i++) {
            cumCC = (int)ceil(r_prob() * max_dist_number);
            if (seq[i] == 'A') {
                auto it = a_qual_dist_first[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'C') {
                auto it = c_qual_dist_first[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'G') {
                auto it = g_qual_dist_first[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'T') {
                auto it = t_qual_dist_first[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else {
                // return random quality less than 10
                read_qual.push_back((int)r_prob() * 10);
            }
        }
        return read_qual;
    }

    vector<int> Empdist::get_read_qual_sep_2(const string& seq) const
    {
        vector<int> read_qual;
        const auto len = static_cast<int>(seq.size());

        if (a_qual_dist_second.size() < len || t_qual_dist_second.size() < len || g_qual_dist_second.size() < len
            || c_qual_dist_second.size() < len) {
            BOOST_LOG_TRIVIAL(fatal) << "Error: The required read length exceeds the "
                                        "length of the read quality profile";
            exit(EXIT_FAILURE);
        }
        int cumCC;
        for (int i = 0; i < len; i++) {
            cumCC = (int)ceil(r_prob() * max_dist_number);
            if (seq[i] == 'A') {
                auto it = a_qual_dist_second[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'C') {
                auto it = c_qual_dist_second[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'G') {
                auto it = g_qual_dist_second[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else if (seq[i] == 'T') {
                auto it = t_qual_dist_second[i].lower_bound(cumCC);
                read_qual.push_back(it->second);
            } else {
                // return random quality less than 10
                read_qual.push_back((int)r_prob() * 10);
            }
        }
        return read_qual;
    }

    void Empdist::read_emp_dist(istream& input, bool is_first)
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
            if (aLine[0] == '#')
                continue;

            if (!_sep_qual && (aLine[0] != '.')) {
                continue;
            }
            if (_sep_qual) {
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
                    BOOST_LOG_TRIVIAL(fatal)
                        << "Fatal error (1): Wrong format of input distribution.";
                    exit(EXIT_FAILURE);
                }
            }

            int t_int;
            vector<int> qual;

            while (ss >> t_int) {
                qual.push_back(t_int);
            }

            getline(input, aLine);
            ss.clear();
            ss.str(aLine);
            ss >> alt_read_pos;
            ss >> read_pos;

            if (read_pos != linenum) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "Fatal error (2): Wrong format of input distribution.";
                exit(EXIT_FAILURE);
            }

            long t_uint;
            vector<long> count;

            while (ss >> t_uint) {
                count.push_back(t_uint);
            }

            if (count.size() != qual.size()) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "Fatal error (3): Wrong format of input distribution.";
                exit(EXIT_FAILURE);
            }

            double denom = static_cast<double>(count[count.size() - 1]) / max_dist_number;
            map<int, int> dist;

            for (size_t i = 0; i < count.size(); i++) {
                auto cc = static_cast<int>(ceil(static_cast<double>(count[i]) / denom));
                dist[cc] = qual[i];
            }
            if (!dist.empty()) {
                linenum++;
                if (!_sep_qual && is_first) {
                    qual_dist_first.push_back(dist);
                } else if (!_sep_qual) {
                    qual_dist_second.push_back(dist);
                } else if (a_flag && is_first) {
                    a_qual_dist_first.push_back(dist);
                } else if (a_flag) {
                    a_qual_dist_second.push_back(dist);
                } else if (t_flag && is_first) {
                    t_qual_dist_first.push_back(dist);
                } else if (t_flag) {
                    t_qual_dist_second.push_back(dist);
                } else if (g_flag && is_first) {
                    g_qual_dist_first.push_back(dist);
                } else if (g_flag) {
                    g_qual_dist_second.push_back(dist);
                } else if (c_flag && is_first) {
                    c_qual_dist_first.push_back(dist);
                } else if (c_flag) {
                    c_qual_dist_second.push_back(dist);
                } else {
                    BOOST_LOG_TRIVIAL(fatal)
                        << "Unexpected Error: Profile was not read in correctly.";
                    exit(EXIT_FAILURE);
                }
            }
        }

        if (_sep_qual) {
            if (a_qual_dist_first.size() != g_qual_dist_first.size() || g_qual_dist_first.size() != c_qual_dist_first.size()
                || c_qual_dist_first.size() != t_qual_dist_first.size()) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "Unexpected Error: Profile was not read in correctly.";
                exit(EXIT_FAILURE);
            }
            if (a_qual_dist_second.size() != g_qual_dist_second.size()
                || g_qual_dist_second.size() != c_qual_dist_second.size()
                || c_qual_dist_second.size() != t_qual_dist_second.size()) {
                BOOST_LOG_TRIVIAL(fatal)
                    << "Unexpected Error: Profile was not read in correctly.";
                exit(EXIT_FAILURE);
            }
        }

        if (linenum == 0) {
            BOOST_LOG_TRIVIAL(fatal) << "Fatal Error: Profile empty!";
            exit(EXIT_FAILURE);
        }
    }

    void Empdist::read_emp_dist(const string& infile, bool is_first)
    {
        ifstream distss(infile.c_str());
        if (!distss) {
            BOOST_LOG_TRIVIAL(fatal)
                << "Fatal Error: Cannot open the distribution file: '" << infile << "'";
            exit(EXIT_FAILURE);
        }
        read_emp_dist(distss, is_first);
    }
    double r_prob()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        return (double)gen() / std::mt19937::max();
    }
}
}