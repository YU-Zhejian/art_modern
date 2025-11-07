#include "art_modern_config.h" // NOLINT: For WITH_BOOST_TIMER, WITH_MPI

#include "art_profile_builder/exe/main_fn.hh"
#include "art_profile_builder/exe/parse_args.hh"
#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include "libam_support/utils/dump_utils.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

using namespace labw::art_modern;

int main(int argc, char** argv)
{
    const auto wall_time_start = std::chrono::high_resolution_clock ::now();
    init_mpi(&argc, &argv);

    if (!is_on_mpi_main_process_or_nompi()) {
        // Only rank 0 does the work
        exit_mpi();
        return EXIT_SUCCESS;
    }
    init_logger();
    init_file_logger("art_profile_builder");
    handle_dumps();

    const APBConfig config = parse_args(argc, argv);

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif

    std::vector<std::thread> threads;

    std::vector<std::shared_ptr<IntermediateEmpDist>> ieds_r1;
    std::vector<std::shared_ptr<IntermediateEmpDist>> ieds_r2;

    for (std::size_t thread_id = 0; thread_id < config.num_threads; ++thread_id) {
        auto this_ied1 = std::make_shared<IntermediateEmpDist>(config.read_length_1);
        ieds_r1.emplace_back(this_ied1);
        auto this_ied2 = std::make_shared<IntermediateEmpDist>(config.read_length_2);
        ieds_r2.emplace_back(this_ied2);
        threads.emplace_back(view_sam, this_ied1, this_ied2, thread_id, config);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    std::size_t total_reads = 0;
    std::size_t total_bases = 0;

    IntermediateEmpDist ied1(config.read_length_1);
    for (const auto& other_ied : ieds_r1) {
        ied1.add(*other_ied);
    }
    ied1.accumulate();
    std::ofstream out1_s { config.output_1_file_path, std::ios::out };
    ied1.write(out1_s, config.is_ob);
    total_reads += ied1.get_total_reads();
    total_bases += ied1.get_total_bases();
    out1_s.close();

    if (config.is_pe) {
        IntermediateEmpDist ied2(config.read_length_2);
        for (const auto& other_ied : ieds_r2) {
            ied2.add(*other_ied);
        }
        ied2.accumulate();
        std::ofstream out2_s { config.output_2_file_path, std::ios::out };
        ied2.write(out2_s, config.is_ob);
        total_reads += ied2.get_total_reads();
        total_bases += ied2.get_total_bases();
        out2_s.close();
    }
    const auto wall_time_end = std::chrono::high_resolution_clock ::now();
    const auto wall_time_diff = wall_time_end - wall_time_start;
    const double wall_time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(wall_time_diff).count();
    BOOST_LOG_TRIVIAL(info) << "Total reads processed: "
                            << to_si(total_reads, 2, static_cast<decltype(total_reads)>(1000))
                            << " speed: " << to_si(static_cast<double>(total_reads) / wall_time_sec, 2, 1000.0)
                            << " reads/s";
    BOOST_LOG_TRIVIAL(info) << "Total bases processed: "
                            << to_si(total_bases, 2, static_cast<decltype(total_bases)>(1000))
                            << " speed: " << to_si(static_cast<double>(total_bases) / wall_time_sec, 2, 1000.0)
                            << " bases/s";

#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time elapsed: " << t.format();
#endif
    BOOST_LOG_TRIVIAL(info) << "Done.";
    exit_mpi();
    return EXIT_SUCCESS;
}
