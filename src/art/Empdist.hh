#pragma once

#include "art_modern_config.h" // NOLINT

#include "art/BuiltinProfile.hh"
#include "art/random_generator.hh"

#include "libam/Dtypes.hh"

#ifdef USE_BTREE_MAP
#include <btree/map.h>
#endif

#include <cstddef>
#include <functional>
#include <istream>
#include <string>
#include <vector>

namespace labw::art_modern {

class Empdist {
public:
#ifdef USE_BTREE_MAP
    using dist_map_type = btree::map<int, int, std::less<>>;
#else
    using dist_map_type = std::map<int, int, std::less<>>;
#endif
    using dist_type = std::vector<dist_map_type>;
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

    Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual, bool is_pe,
        std::size_t read_len);

    Empdist(const BuiltinProfile& builtin_profile, bool sep_qual, bool is_pe, std::size_t read_len);
    void get_read_qual(std::vector<am_qual_t>& qual, int len, Rprob& rprob, bool first = true) const;
    void get_read_qual_sep_1(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    void get_read_qual_sep_2(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    void shift_all_emp(bool sep_flag, int q_shift_1, int q_shift_2, int min_qual, int max_qual);

private:
    void read_emp_dist_(const std::string& infile, bool is_first);
    void read_emp_dist_(std::istream& input, bool is_first);
    void validate_() const;
    void print_() const;
    bool sep_qual_;
    bool is_pe_;
    std::size_t read_len_;
};

} // namespace labw::art_modern
