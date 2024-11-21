#include "art_modern_config.h"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <fstream>
#include <iostream>
#include <thread>

#include "ArtConstants.hh"
#include "ArtJobPool.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/FastaStreamBatcher.hh"
#include "fasta/Pbsim3TranscriptBatcher.hh"
#include "global_variables.hh"
#include "main_fn.hh"

namespace labw::art_modern {

void print_banner()
{
    BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina (" << ART_PROGRAM_NAME << ")";
    BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
    BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
    BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejianyu@intl.zju.edu.cn>";
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
                    SimulationJob sj(new InMemoryFastaFetch(fa_view), coverage_info, job_id);
                    ArtJobExecutor aje(std::move(sj), art_params);
                    job_pool.add(std::move(aje));
                }
            } else {
                std::ifstream fasta_stream(art_params.input_file_name);
                FastaStreamBatcher fsb(art_params.batch_size, fasta_stream);
                while (true) {
                    auto fa_view = fsb.fetch();

                    if (fa_view.num_seqs() == 0) {
                        break;
                    }
                    thread_local SimulationJob sj(new InMemoryFastaFetch(fa_view), coverage_info, ++job_id);
                    ArtJobExecutor aje(std::move(sj), art_params);
                    job_pool.add(std::move(aje));
                }
                fasta_stream.close();
            }
        }
        if (art_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TEMPLATE) {
            if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                std::ifstream pbsim3_template_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_params.batch_size, pbsim3_template_stream);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.first.num_seqs() == 0) {
                        break;
                    }
                    thread_local SimulationJob sj(new InMemoryFastaFetch(fa_view.first), fa_view.second, ++job_id);
                    ArtJobExecutor aje(std::move(sj), art_params);
                    job_pool.add(std::move(aje));
                }
                pbsim3_template_stream.close();
            } else {
                std::ifstream pbsim3_template_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_params.batch_size, pbsim3_template_stream);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.first.num_seqs() == 0) {
                        break;
                    }
                    thread_local SimulationJob sj(new InMemoryFastaFetch(fa_view.first), fa_view.second, ++job_id);
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
} // namespace labw::art_modern
