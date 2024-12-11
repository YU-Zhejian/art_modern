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
    using dist_type = std::vector<std::map<int, int, std::less<>>>;
    dist_type qual_dist_first;
    dist_type qual_dist_second;

    dist_type a_qual_dist_first;
    dist_type t_qual_dist_first;
    dist_type g_qual_dist_first;
    dist_type c_qual_dist_first;

    dist_type a_qual_dist_second;
    dist_type t_qual_dist_second;
    dist_type g_qual_dist_second;
    dist_type c_qual_dist_second;

    Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual);
    void get_read_qual(std::vector<int>& qual, int len, Rprob& rprob, bool first = true) const;
    void get_read_qual_sep_1(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const;
    void get_read_qual_sep_2(std::vector<int>& qual, const std::string& seq, Rprob& rprob) const;
    void shift_all_emp(
        const bool sep_flag, const int q_shift_1, const int q_shift_2, const int min_qual, const int max_qual);

private:
    void read_emp_dist_(const std::string& infile, bool is_first);
    void read_emp_dist_(std::istream& input, bool is_first);
    bool sep_qual_;
};

} // namespace labw::art_modern // namespace labw