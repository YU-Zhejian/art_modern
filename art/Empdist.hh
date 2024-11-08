#pragma once
#include <array>
#include <istream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "random_generator.hh"

namespace labw::art_modern {

class Empdist {
public:
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
    void get_read_qual(std::vector<int>& qual, int len, Rprob& rprob, bool first = true) const;
    void get_read_qual_sep_1(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const;
    void get_read_qual_sep_2(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const;

private:
    void read_emp_dist_(const std::string& infile, bool is_first);
    void read_emp_dist_(std::istream& input, bool is_first);
    bool sep_qual_;
};

} // namespace labw::art_modern // namespace labw