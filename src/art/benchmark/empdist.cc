#include <chrono>

#include "art/lib/ArtConstants.hh"
#include "art/lib/Empdist.hh"
#include "art/lib/Rprob.hh"

#include "benchmark/benchmark_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"

#include <boost/log/trivial.hpp>

#include <chrono>
#include <cstdlib>
#include <vector>

using namespace labw::art_modern;
constexpr std::size_t N_REPLICA = 500UL;
constexpr std::size_t N_TIMES = M_SIZE * static_cast<std::size_t>(10);

int main() {
    Empdist const empdist{DEFAULT_ERR_PROFILE, false, true, false };

    Rprob rp(200.0, 10.0, 10, 150);
    std::vector<am_qual_t>qual (150) ;
    std::vector<std::size_t> times;

    for (std::size_t i = 0; i < N_REPLICA; i++)
    {
        const auto start = std::chrono::high_resolution_clock::now();
                for (std::size_t j = 0; j < N_TIMES; j++){
            empdist.get_read_qual( qual, rp, true );
            empdist.get_read_qual( qual, rp, false);
        }
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = end - start;
        times.emplace_back(std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());
    }
    BOOST_LOG_TRIVIAL(info) << "Empdist::get_read_qual  : " << describe( times, true);
    // Intel:
    // [2025-12-18 23:15:30.327662] [0x00007e72041eb200] [info]    Empdist::get_read_qual  : gmean: 35,650,559; mean/sd: 35,743,146/2,945,319
    // GCC:
    // [2025-12-18 23:16:54.303635] [0x0000712f376021c0] [info]    Empdist::get_read_qual  : gmean: 57,546,368; mean/sd: 59,443,214/20,806,297
    return EXIT_SUCCESS;
}
