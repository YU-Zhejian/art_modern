// TODO:
// 1. Add argument parsing

#include "art_profile_builder/exe/APBConfig.hh"
#include "art_profile_builder/exe/IntermediateEmpDist.hh"
#include "art_profile_builder/exe/main_fn.hh"
#include "art_profile_builder/exe/parse_args.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "libam_support/utils/mpi_utils.hh"

using namespace labw::art_modern;

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);

#ifdef WITH_MPI
    if (mpi_rank() != MPI_MAIN_RANK_STR) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
#endif

#ifdef WITH_BOOST_TIMER
    boost::timer::cpu_timer timer;
#endif

    const APBConfig config = parse_args(argc, argv);

    std::vector<std::thread> threads;
    threads.reserve(config.num_threads);

    std::vector<std::shared_ptr<IntermediateEmpDist>> ieds_r1;
    std::vector<std::shared_ptr<IntermediateEmpDist>> ieds_r2;

    for (std::size_t thread_id = 0; thread_id < config.num_threads; ++thread_id) {
        auto this_ied1 = std::make_shared<IntermediateEmpDist>(config.read_length);
        ieds_r1.emplace_back(this_ied1);

        auto this_ied2 = std::make_shared<IntermediateEmpDist>(config.read_length);
        ieds_r2.emplace_back(this_ied2);

        threads.emplace_back(view_sam, this_ied1, this_ied2, thread_id, config);
    }
    for (auto& t : threads) {
        t.join();
    }

    IntermediateEmpDist ied1(config.read_length);
    for (const auto& other_ied : ieds_r1) {
        ied1.add(*other_ied);
    }
    ied1.accumulate();
    ied1.write(std::cout);

    if (config.is_pe) {
        IntermediateEmpDist ied2(config.read_length);
        for (const auto& other_ied : ieds_r2) {
            ied2.add(*other_ied);
        }
        ied2.accumulate();
        ied2.write(std::cout);
    }

    BOOST_LOG_TRIVIAL(info) << "Done.";

#ifdef WITH_BOOST_TIMER
    timer.stop();
    BOOST_LOG_TRIVIAL(info) << "Time elapsed: " << timer.format();
#endif
    // FIXME: Redo this analysis
    // Threads: 1: Time elapsed:  11.390000s wall, 11.170000s user + 0.210000s system = 11.380000s CPU (99.9%)
    // Threads: 2: Time elapsed:  7.330000s wall, 14.270000s user + 0.320000s system = 14.590000s CPU (199.0%)

    exit_mpi();
    return EXIT_SUCCESS;
}
