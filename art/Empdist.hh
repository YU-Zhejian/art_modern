#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "ArtConstants.hh"

namespace labw {
namespace art_modern {

    double r_prob();

    class Empdist {
    public:
        static const int max_dist_number = 1000000;
        std::vector<std::map<int, int>> qual_dist_first;
        std::vector<std::map<int, int>> qual_dist_second;

        std::vector<std::map<int, int>> a_qual_dist_first;
        std::vector<std::map<int, int>> t_qual_dist_first;
        std::vector<std::map<int, int>> g_qual_dist_first;
        std::vector<std::map<int, int>> c_qual_dist_first;

        std::vector<std::map<int, int>> a_qual_dist_second;
        std::vector<std::map<int, int>> t_qual_dist_second;
        std::vector<std::map<int, int>> g_qual_dist_second;
        std::vector<std::map<int, int>> c_qual_dist_second;

        Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual);
        std::vector<int> get_read_qual(int len, bool first = true) const;
        std::vector<int> get_read_qual_sep_1(const std::string& seq) const;
        std::vector<int> get_read_qual_sep_2(const std::string& seq) const;

    private:
        bool _sep_qual;

        void read_emp_dist(const std::string& infile, bool is_first);
        void read_emp_dist(std::istream& input, bool is_first);
    };

} // namespace art_modern
} // namespace labw