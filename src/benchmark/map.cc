/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/ds/GslDiscreteDistribution.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/accumulators/accumulators_fwd.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/map.hpp>
#include <boost/log/trivial.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace labw::art_modern;

namespace {
/**
 * From base 2 of MiSeqv3L250R1.
 */
const std::vector<am_qual_t> empdist_quals { 2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38 };
const std::vector<am_qual_count_t> empdist_counts { 1'005'307, 1'005'331, 1'052'941, 1'212'250, 1'280'823, 1'319'284,
    1'319'306, 1'322'792, 1'348'204, 1'378'493, 1'402'848, 1'437'193, 1'470'236, 1'493'447, 1'558'624, 1'588'409,
    1'687'992, 1'729'588, 1'809'750, 1'898'693, 1'926'134, 2'010'366, 2'050'530, 2'111'370, 2'160'093, 2'209'413,
    2'245'481, 2'271'272, 2'618'106, 2'797'257, 2'970'924, 3'577'430, 6'455'505 };
const auto DIST_END = empdist_counts.back();
constexpr am_qual_count_t MAX_DIST_NUMBER = 1'000'000;

constexpr am_readnum_t NUM_TRIALS = K_SIZE << 8;
constexpr am_readnum_t READ_LEN = K_SIZE;
std::mt19937 gen { 0 };
std::uniform_int_distribution<am_qual_count_t> one_to_dist_end(1, DIST_END);
std::uniform_int_distribution<am_qual_count_t> one_to_max_dist_number(1, MAX_DIST_NUMBER);
} // namespace

class SlimEmpDist {
public:
    virtual void gen_qualities(std::vector<am_qual_t>& qual) = 0;
    SlimEmpDist() = default;
    DELETE_MOVE(SlimEmpDist)
    DELETE_COPY(SlimEmpDist)
    virtual ~SlimEmpDist() = default;
};

class SlimEmpDistNop : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistNop)
    DELETE_COPY(SlimEmpDistNop)
    SlimEmpDistNop() = default;
    ~SlimEmpDistNop() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (int i = 0; i < READ_LEN; i++) {
            gen();
        }
    }
};
class SlimEmpDistOld : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistOld)
    DELETE_COPY(SlimEmpDistOld)
    SlimEmpDistOld(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
    {
        const auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;
        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
        }
    }
    ~SlimEmpDistOld() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual) override
    {
        std::vector<am_qual_count_t> count;
        count.resize(READ_LEN);
        std::generate(count.begin(), count.end(), [] { return one_to_max_dist_number(gen); });
        for (std::size_t i = 0; i < qual.size(); i++) {
            qual[i] = dist_.lower_bound(count[i])->second;
        }
    }

private:
    std::map<am_qual_count_t, am_qual_t, std::less<>> dist_;
};
class SlimEmpDistStdDiscrete : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistStdDiscrete)
    DELETE_COPY(SlimEmpDistStdDiscrete)
    SlimEmpDistStdDiscrete(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
        : qual_(qual)
    {
        std::vector<double> init_list;
        double prev = 0;
        for (const auto i : count) {
            init_list.emplace_back(i - prev);
            prev = i;
        }
        rd_ = std::discrete_distribution<std::size_t> { init_list.begin(), init_list.end() };
    }
    ~SlimEmpDistStdDiscrete() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (auto& i : qual) {
            i = qual_[rd_(gen)];
        }
    }

private:
    std::vector<am_qual_t> qual_;
    std::discrete_distribution<std::size_t> rd_;
};
class SlimEmpDistBoostDiscrete : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistBoostDiscrete)
    DELETE_COPY(SlimEmpDistBoostDiscrete)
    SlimEmpDistBoostDiscrete(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
        : qual_(qual)
    {
        std::vector<double> init_list;
        double prev = 0;
        for (const int i : count) {
            init_list.emplace_back(i - prev);
            prev = i;
        }
        rd_ = boost::random::discrete_distribution<std::size_t> { init_list.begin(), init_list.end() };
    }
    ~SlimEmpDistBoostDiscrete() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (auto& i : qual) {
            i = qual_[rd_(gen)];
        }
    }

private:
    std::vector<am_qual_t> qual_;
    boost::random::discrete_distribution<std::size_t> rd_;
};
class SlimEmpDistGslDiscrete : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistGslDiscrete)
    DELETE_COPY(SlimEmpDistGslDiscrete)
    SlimEmpDistGslDiscrete(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
        : qual_(qual)
    {
        std::vector<double> init_list;
        double prev = 0;
        for (const int i : count) {
            init_list.emplace_back((i - prev));
            prev = i;
        }
        rd_ = GslDiscreteDistribution<double>(init_list);
    }
    ~SlimEmpDistGslDiscrete() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (auto& i : qual) {
            i = qual_[rd_(rnd_01_(gen))];
        }
    }

private:
    std::vector<am_qual_t> qual_;
    GslDiscreteDistribution<double> rd_;
    boost::uniform_01<double> rnd_01_;
};
class SlimEmpDistGslDiscreteInt : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistGslDiscreteInt)
    DELETE_COPY(SlimEmpDistGslDiscreteInt)
    SlimEmpDistGslDiscreteInt(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
        : qual_(qual)
        , rnd_(0, count.size() - 1)
    {
        std::vector<int64_t> init_list;
        int64_t prev = 0;
        for (const int i : count) {
            init_list.emplace_back((i - prev));
            prev = i;
        }
        rd_ = GslDiscreteIntDistribution<int64_t>(init_list);
    }
    ~SlimEmpDistGslDiscreteInt() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (auto& i : qual) {
            i = qual_[rd_(rnd_(gen))];
        }
    }

private:
    std::vector<am_qual_t> qual_;
    GslDiscreteIntDistribution<int64_t> rd_;
    boost::random::uniform_int_distribution<int64_t> rnd_;
};
class SlimEmpDistGslDiscreteInterpolated : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistGslDiscreteInterpolated)
    DELETE_COPY(SlimEmpDistGslDiscreteInterpolated)
    SlimEmpDistGslDiscreteInterpolated(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
        : qual_offset_(qual.front())
    {
        am_qual_t prev_qual = qual_offset_;
        std::vector<double> init_list;
        double prev = 0;
        for (std::size_t i = 0; i < count.size(); i++) {
            if (i != 0) {
                prev_qual++;
            }
            while (prev_qual != qual[i]) {
                init_list.emplace_back(0);
                prev_qual++;
            }
            init_list.emplace_back(count[i] - prev);
            prev = count[i];
        }
        rd_ = GslDiscreteDistribution<double>(init_list);
    }
    ~SlimEmpDistGslDiscreteInterpolated() override = default;

    void gen_qualities([[maybe_unused]] std::vector<am_qual_t>& qual) override
    {
        for (auto& i : qual) {
            i = rd_(rnd_01_(gen)) + qual_offset_;
        }
    }

private:
    am_qual_t qual_offset_;
    GslDiscreteDistribution<double> rd_;
    boost::uniform_01<double> rnd_01_;
};
class SlimEmpDistUsingBoostMap : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistUsingBoostMap)
    DELETE_COPY(SlimEmpDistUsingBoostMap)
    SlimEmpDistUsingBoostMap(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
    {
        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[count[i]] = qual[i];
        }
    }
    ~SlimEmpDistUsingBoostMap() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual) override
    {
        std::vector<am_qual_count_t> count;
        count.resize(READ_LEN);
        std::generate(count.begin(), count.end(), [] { return one_to_dist_end(gen); });
        for (std::size_t i = 0; i < qual.size(); i++) {
            qual[i] = dist_.lower_bound(count[i])->second;
        }
    }

private:
    boost::container::map<am_qual_count_t, am_qual_t, std::less<>> dist_;
};

class SlimEmpDistUsingBoostFlatMap : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistUsingBoostFlatMap)
    DELETE_COPY(SlimEmpDistUsingBoostFlatMap)
    SlimEmpDistUsingBoostFlatMap(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
    {
        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[count[i]] = qual[i];
        }
    }
    ~SlimEmpDistUsingBoostFlatMap() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual) override
    {
        std::vector<am_qual_count_t> count;
        count.resize(READ_LEN);
        std::generate(count.begin(), count.end(), [] { return one_to_dist_end(gen); });
        for (std::size_t i = 0; i < qual.size(); i++) {
            qual[i] = dist_.lower_bound(count[i])->second;
        }
    }

private:
    boost::container::flat_map<am_qual_count_t, am_qual_t, std::less<>> dist_;
};

class SlimEmpDistUsingStdMap : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistUsingStdMap)
    DELETE_COPY(SlimEmpDistUsingStdMap)
    SlimEmpDistUsingStdMap(const std::vector<am_qual_t>& qual, const std::vector<am_qual_count_t>& count)
    {
        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[count[i]] = qual[i];
        }
    }
    ~SlimEmpDistUsingStdMap() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual) override
    {
        std::vector<am_qual_count_t> count;
        count.resize(READ_LEN);
        std::generate(count.begin(), count.end(), [] { return one_to_dist_end(gen); });
        for (std::size_t i = 0; i < qual.size(); i++) {
            qual[i] = dist_.lower_bound(count[i])->second;
        }
    }

private:
    std::map<am_qual_count_t, am_qual_t, std::less<>> dist_;
};

namespace {
void bench(std::unique_ptr<SlimEmpDist> empdist, const std::string& name)
{
    boost::accumulators::accumulator_set<double,
        boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance>>
        acc;

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::vector<am_qual_t> qual;
    qual.resize(READ_LEN);
    start = std::chrono::high_resolution_clock::now();
    for (am_readnum_t i = 0; i < NUM_TRIALS; i++) {
        empdist->gen_qualities(qual);
        for (const auto value : qual) {
            acc(value);
        }
    }
    end = std::chrono::high_resolution_clock::now();

    // Extract the mean and standard deviation
    const double mean = boost::accumulators::mean(acc);
    const double sd = std::sqrt(boost::accumulators::variance(acc));
    BOOST_LOG_TRIVIAL(info) << name << ": "
                            << format_with_commas(
                                   std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count())
                            << " ns. Q: mean=" << mean << " sd=" << sd;
}
} // namespace

int main()
{
    bench(std::make_unique<SlimEmpDistNop>(), "nop");
    bench(std::make_unique<SlimEmpDistOld>(empdist_quals, empdist_counts), "old");
    bench(std::make_unique<SlimEmpDistUsingStdMap>(empdist_quals, empdist_counts), "std::map");
    bench(std::make_unique<SlimEmpDistUsingBoostMap>(empdist_quals, empdist_counts), "boost::map");
    bench(std::make_unique<SlimEmpDistUsingBoostFlatMap>(empdist_quals, empdist_counts), "boost::flat_map");
    bench(std::make_unique<SlimEmpDistStdDiscrete>(empdist_quals, empdist_counts), "std::dd");
    bench(std::make_unique<SlimEmpDistBoostDiscrete>(empdist_quals, empdist_counts), "boost::dd");
    bench(std::make_unique<SlimEmpDistGslDiscrete>(empdist_quals, empdist_counts), "GSL::dd (float)");
    // Not working since accuracy lost?
    bench(std::make_unique<SlimEmpDistGslDiscreteInt>(empdist_quals, empdist_counts), "GSL::dd (int)");
    bench(std::make_unique<SlimEmpDistGslDiscreteInterpolated>(empdist_quals, empdist_counts), "GSL::dd (+)");
    return EXIT_SUCCESS;
}
