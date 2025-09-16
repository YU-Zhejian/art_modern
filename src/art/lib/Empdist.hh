#pragma once

#include "art_modern_config.h" // NOLINT

#include "art/lib/BuiltinProfile.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/class_macros_utils.hh"

#if defined(USE_WALKER_QUALGEN)
#include "libam_support/ds/GslDiscreteDistribution.hh"
#endif

#include <cstddef>
#include <functional>
#include <istream>
#include <map>
#include <string>
#include <vector>

namespace labw::art_modern {

class Empdist {
public:
    using dist_map_type = std::map<am_qual_count_t, am_qual_t, std::less<>>;
    using dist_type = std::vector<dist_map_type>;

private:
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
#ifdef USE_WALKER_QUALGEN
    class SlimEmpDistGslDiscrete {
    public:
        DEFAULT_COPY(SlimEmpDistGslDiscrete)
        DEFAULT_MOVE(SlimEmpDistGslDiscrete)
        explicit SlimEmpDistGslDiscrete(const dist_map_type& dist)
        {
            std::vector<double> count;
            for (const auto& [this_count, this_qual] : dist) {
                qual_.emplace_back(this_qual);
                count.emplace_back(static_cast<double >(this_count));
            }
            std::vector<double> init_list;
            double prev = 0;
            for (const auto i : count) {
                init_list.emplace_back(i - prev);
                prev = i;
            }
            rd_ = GslDiscreteDistribution<double>(init_list);
        }
        ~SlimEmpDistGslDiscrete() = default;

        [[nodiscard]] am_qual_t gen_qual(const double u) const { return qual_[rd_(u)]; }

    private:
        std::vector<am_qual_t> qual_;
        GslDiscreteDistribution<double> rd_;
    };

    using dist_idx_type = std::vector<SlimEmpDistGslDiscrete>;
    dist_idx_type qual_dist_first_idx;
    dist_idx_type qual_dist_second_idx;
    dist_idx_type a_qual_dist_first_idx;
    dist_idx_type t_qual_dist_first_idx;
    dist_idx_type g_qual_dist_first_idx;
    dist_idx_type c_qual_dist_first_idx;
    dist_idx_type a_qual_dist_second_idx;
    dist_idx_type t_qual_dist_second_idx;
    dist_idx_type g_qual_dist_second_idx;
    dist_idx_type c_qual_dist_second_idx;
#endif

public:
    Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual, bool is_pe,
        std::size_t read_len);

    Empdist(const BuiltinProfile& builtin_profile, bool sep_qual, bool is_pe, std::size_t read_len);
    void get_read_qual(std::vector<am_qual_t>& qual, Rprob& rprob, bool first = true) const;
    void get_read_qual_sep_1(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    void get_read_qual_sep_2(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    void shift_all_emp(bool sep_flag, am_qual_t q_shift_1, am_qual_t q_shift_2, am_qual_t min_qual, am_qual_t max_qual);
    void index();
    void log() const;

private:
    void read_emp_dist_(const std::string& infile, bool is_first);
    void read_emp_dist_(std::istream& input, bool is_first);
    void validate_() const;
    const bool sep_qual_;
    const bool is_pe_;
    const std::size_t read_len_;
};

} // namespace labw::art_modern
