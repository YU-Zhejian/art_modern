#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <btree/map.h>

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
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
const std::vector<am_qual_dist_t> empdist_counts { 1005307, 1005331, 1052941, 1212250, 1280823, 1319284, 1319306,
    1322792, 1348204, 1378493, 1402848, 1437193, 1470236, 1493447, 1558624, 1588409, 1687992, 1729588, 1809750, 1898693,
    1926134, 2010366, 2050530, 2111370, 2160093, 2209413, 2245481, 2271272, 2618106, 2797257, 2970924, 3577430,
    6455505 };
constexpr am_qual_dist_t MAX_DIST_NUMBER = 1000000;
constexpr am_readnum_t NUM_TRIALS = K_SIZE << 8;
constexpr am_readnum_t READ_LEN = K_SIZE;

} // namespace

class SlimEmpDist {
public:
    virtual void gen_qualities(std::vector<am_qual_t>& qual, const std::vector<am_qual_dist_t>& count) const = 0;
    SlimEmpDist() = default;
    DELETE_MOVE(SlimEmpDist)
    DELETE_COPY(SlimEmpDist)
    virtual ~SlimEmpDist() = default;
};

class SlimEmpDistUsingStdMap : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistUsingStdMap)
    DELETE_COPY(SlimEmpDistUsingStdMap)
    SlimEmpDistUsingStdMap(const std::vector<am_qual_t>& qual, const std::vector<am_qual_dist_t>& count)
    {
        const auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;

        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
        }
    }
    ~SlimEmpDistUsingStdMap() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual, const std::vector<am_qual_dist_t>& count) const override
    {
        for (std::size_t i = 0; i < count.size(); i++) {
            qual[i] = dist_.lower_bound(qual[i])->second;
        }
    }

private:
    std::map<am_qual_t, am_qual_dist_t, std::less<>> dist_;
};

class SlimEmpDistUsingBTreeMap : public SlimEmpDist {
public:
    DELETE_MOVE(SlimEmpDistUsingBTreeMap)
    DELETE_COPY(SlimEmpDistUsingBTreeMap)
    SlimEmpDistUsingBTreeMap(const std::vector<am_qual_t>& qual, const std::vector<am_qual_dist_t>& count)
    {
        const auto denom = static_cast<double>(count.back()) / MAX_DIST_NUMBER;

        for (decltype(count.size()) i = 0; i < count.size(); i++) {
            dist_[static_cast<int>(std::ceil(static_cast<double>(count[i]) / denom))] = qual[i];
        }
    }
    ~SlimEmpDistUsingBTreeMap() override = default;

    void gen_qualities(std::vector<am_qual_t>& qual, const std::vector<am_qual_dist_t>& count) const override
    {
        for (std::size_t i = 0; i < count.size(); i++) {
            qual[i] = dist_.lower_bound(qual[i])->second;
        }
    }

private:
    btree::map<am_qual_t, am_qual_dist_t, std::less<>> dist_;
};

namespace {
void bench(std::unique_ptr<SlimEmpDist> empdist, const std::string& name, const std::vector<am_qual_dist_t>& count)
{
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::vector<am_qual_t> qual;
    qual.resize(count.size());
    start = std::chrono::high_resolution_clock::now();
    for (am_readnum_t i = 0; i < NUM_TRIALS; i++) {
        empdist->gen_qualities(qual, count);
    }
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << name << ": " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
                            << " ns";
}
} // namespace

int main()
{
    std::vector<am_qual_dist_t> count;
    count.resize(READ_LEN);
    std::mt19937 gen { std::random_device {}() };
    std::generate(count.begin(), count.end(),
        std::bind(std::uniform_int_distribution<am_qual_dist_t>(1, MAX_DIST_NUMBER), std::ref(gen)));

    bench(std::make_unique<SlimEmpDistUsingStdMap>(empdist_quals, empdist_counts), "std::map", empdist_counts);
    bench(std::make_unique<SlimEmpDistUsingBTreeMap>(empdist_quals, empdist_counts), "btree::map", empdist_counts);
    return EXIT_SUCCESS;
}