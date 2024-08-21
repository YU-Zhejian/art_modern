#include "ArtJobPool.hh"
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>

#include "ArtConstants.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/FastaStreamBatcher.hh"
#include "fasta/Pbsim3TranscriptBatcher.hh"
#include "main_fn.hh"
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/stacktrace.hpp>
#include <csignal>
#include <iostream>
#define DUMP_FILENAME "./backtrace.dump"

using namespace std;
namespace logging = boost::log;

namespace labw {
namespace art_modern {

    void init_logger()
    {
#ifndef CEU_CM_IS_DEBUG
        logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
#endif
    }

    void my_signal_handler(int signum)
    {
        ::signal(signum, SIG_DFL);
        boost::stacktrace::safe_dump_to(DUMP_FILENAME);
        ::raise(SIGABRT);
    }

    void handle_dumps()
    {
        ::signal(SIGSEGV, &my_signal_handler);
        ::signal(SIGABRT, &my_signal_handler);
        if (boost::filesystem::exists(DUMP_FILENAME)) {
            // there is a backtrace
            std::ifstream ifs(DUMP_FILENAME);

            boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
            std::cout << "Previous run crashed:\n" << st << std::endl;

            // cleaning up
            ifs.close();
            boost::filesystem::remove(DUMP_FILENAME);
        }
    }

    void print_banner()
    {
        BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina (" << ART_PROGRAM_NAME << ")";
        BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
        BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
        BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejian.23@intl.zju.edu.cn>";
#ifdef CEU_CM_IS_DEBUG
        BOOST_LOG_TRIVIAL(info) << "Debugging functions enabled.";
#endif
    }

    void generate_all(const ArtParams& art_params)
    {
        int job_id = 0;
        ArtJobPool job_pool(art_params);
        if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
            // Coverage-based parallelism
            auto coverage_info = art_params.coverage_info.div(art_params.parallel);
            for (int i = 0; i < art_params.parallel; ++i) {
                BaseFastaFetch* fetch;
                if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                    fetch = new InMemoryFastaFetch(art_params.input_file_name);
                } else {
                    fetch = new FaidxFetch(art_params.input_file_name);
                }
                SimulationJob sj(fetch, coverage_info, ++job_id);
                ArtJobExecutor aje(std::move(sj), art_params);
                job_pool.add(std::move(aje));
            }
        } else {
            // Batch-based parallelism
            if (art_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {
                auto const& coverage_info = art_params.coverage_info;
                if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                    InMemoryFastaStreamBatcher fsb(
                        static_cast<int>(art_params.fasta_fetch->num_seqs() / art_params.parallel + 1),
                        art_params.fasta_fetch);
                    while (true) {
                        auto fa_view = fsb.fetch();

                        if (fa_view.num_seqs() == 0) {
                            break;
                        }
                        job_id += 1;
                        SimulationJob sj(new InMemoryFastaFetch(fa_view), coverage_info,  job_id);
                        ArtJobExecutor aje(std::move(sj), art_params);
                        std::cout << "Simulation job " << job_id << endl;
                        job_pool.add(std::move(aje));
                    }
                } else {
                    std::ifstream fasta_stream(art_params.input_file_name);
                    FastaStreamBatcher fsb(1000, fasta_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();

                        if (fa_view.num_seqs() == 0) {
                            break;
                        }
                        thread_local SimulationJob sj (new InMemoryFastaFetch(fa_view), coverage_info, ++job_id);
                        ArtJobExecutor aje(std::move(sj), art_params);
                        job_pool.add(std::move(aje));
                    }
                    fasta_stream.close();
                }
            }
            if (art_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
                if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                    std::ifstream pbsim3_template_stream(art_params.input_file_name);
                    Pbsim3TranscriptBatcher fsb(1000, pbsim3_template_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();
                        if (fa_view.first.num_seqs() == 0) {
                            break;
                        }
                        thread_local SimulationJob sj(  new InMemoryFastaFetch(fa_view.first), fa_view.second, ++job_id );
                        ArtJobExecutor aje(std::move(sj), art_params);
                        job_pool.add(std::move(aje));
                    }
                    pbsim3_template_stream.close();
                } else {
                    std::ifstream pbsim3_template_stream(art_params.input_file_name);
                    Pbsim3TranscriptBatcher fsb(1000, pbsim3_template_stream);
                    while (true) {
                        auto fa_view = fsb.fetch();
                        if (fa_view.first.num_seqs() == 0) {
                            break;
                        }
                        thread_local SimulationJob sj( new InMemoryFastaFetch(fa_view.first), fa_view.second, ++job_id );
                        ArtJobExecutor aje(std::move(sj), art_params);
                        job_pool.add(std::move(aje));
                    }
                    pbsim3_template_stream.close();
                }
            }
        }
        job_pool.stop();
        art_params.out_dispatcher->close();
        delete art_params.out_dispatcher;
        delete art_params.fasta_fetch;
    }
} // namespace art_modern
} // namespace labw
