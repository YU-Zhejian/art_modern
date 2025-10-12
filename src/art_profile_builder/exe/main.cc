#include "art_profile_builder/exe/main_fn.hh"
#include "art_profile_builder/exe/parse_args.hh"
#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

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
    const APBConfig config = parse_args(argc, argv);

#ifdef WITH_BOOST_TIMER
    boost::timer::cpu_timer timer;
#endif

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
    std::ofstream out1_s { config.output_1_file_path, std::ios::out };
    ied1.write(out1_s, config.is_ob);
    out1_s.close();

    if (config.is_pe) {
        IntermediateEmpDist ied2(config.read_length);
        for (const auto& other_ied : ieds_r2) {
            ied2.add(*other_ied);
        }
        ied2.accumulate();
        std::ofstream out2_s { config.output_2_file_path, std::ios::out };
        ied2.write(out2_s, config.is_ob);
        out2_s.close();
    }

    BOOST_LOG_TRIVIAL(info) << "Done.";

#ifdef WITH_BOOST_TIMER
    timer.stop();
    BOOST_LOG_TRIVIAL(info) << "Time elapsed: " << timer.format();
#endif
    exit_mpi();
    return EXIT_SUCCESS;
}
