#pragma once
#include <istream>
#include <map>
#include <string>
#include <vector>

#include "random_generator.hh"

namespace labw {
namespace art_modern {

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

        /**
         * FIXME: TO BE REMOVED
         */
        Empdist();
        Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual);
        std::vector<int> get_read_qual(int len, Rprob& rprob, bool first = true) const;
        std::vector<int> get_read_qual_sep_1(const std::string& seq, Rprob& rprob) const;
        std::vector<int> get_read_qual_sep_2(const std::string& seq, Rprob& rprob) const;

    private:
        bool _sep_qual;

        void read_emp_dist(const std::string& infile, bool is_first);
        void read_emp_dist(std::istream& input, bool is_first);
    };

} // namespace art_modern
} // namespace labw