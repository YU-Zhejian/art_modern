#include "art_modern_config.h"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <fstream>
#include <iostream>

#include "ArtConstants.hh"
#include "ArtJobPool.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/FastaStreamBatcher.hh"
#include "fasta/Pbsim3TranscriptBatcher.hh"
#include "main_fn.hh"
#include "out/OutputDispatcher.hh"

namespace labw::art_modern {

void print_banner()
{
    BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina (" << ART_PROGRAM_NAME << " v. " ART_MODERN_VERSION << ")";
    BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
    BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
    BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejianyu@intl.zju.edu.cn>";
#ifdef CEU_CM_IS_DEBUG
    BOOST_LOG_TRIVIAL(info) << "Debugging functions enabled.";
#endif
}

void generate_wgs(const ArtParams& art_params)
{
    const auto out_dispatcher_factory = get_output_dispatcher_factory();
    BaseReadOutput* out_dispatcher;
    int job_id = 0;
    ArtJobPool job_pool(art_params);

    // Coverage-based parallelism
    const auto& coverage_info = art_params.coverage_info.div(art_params.parallel);

    BaseFastaFetch* fetch;
    if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
        fetch = new InMemoryFastaFetch(art_params.input_file_name);
    } else {
        fetch = new FaidxFetch(art_params.input_file_name);
    }
    out_dispatcher = out_dispatcher_factory.create(art_params.vm, fetch, art_params.args);
    if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
        for (int i = 0; i < art_params.parallel; ++i) {
            SimulationJob sj(fetch, coverage_info, ++job_id, false);
            ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
            job_pool.add(std::move(aje));
        }
    } else {
        for (int i = 0; i < art_params.parallel; ++i) {
            BaseFastaFetch* thread_fetch = new InMemoryFastaFetch(art_params.input_file_name);
            SimulationJob sj(thread_fetch, coverage_info, ++job_id, true);
            ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
            job_pool.add(std::move(aje));
        }
    }
    job_pool.stop();
    BOOST_LOG_TRIVIAL(info) << "Job pool stopped";
    delete fetch;
    out_dispatcher->close();
    delete out_dispatcher;
}

void generate_all(const ArtParams& art_params)
{
    if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
        generate_wgs(art_params);
    } else {
        const auto out_dispatcher_factory = get_output_dispatcher_factory();
        BaseReadOutput* out_dispatcher;
        int job_id = 0;
        ArtJobPool job_pool(art_params);
        // Batch-based parallelism
        if (art_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {
            auto const& coverage_info = art_params.coverage_info;
            if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                InMemoryFastaFetch fetch(art_params.input_file_name);
                InMemoryFastaStreamBatcher fsb(static_cast<int>(fetch.num_seqs() / art_params.parallel + 1), &fetch);
                out_dispatcher = out_dispatcher_factory.create(art_params.vm, &fetch, art_params.args);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.num_seqs() == 0) {
                        break;
                    }
                    job_id += 1;
                    SimulationJob sj(new InMemoryFastaFetch(std::move(fa_view)), coverage_info, job_id, true);
                    ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
                    job_pool.add(std::move(aje));
                }
                // FIXME: fetch not closed!
            } else {
                std::ifstream fasta_stream(art_params.input_file_name);
                FastaStreamBatcher fsb(art_params.batch_size, fasta_stream);
                out_dispatcher
                    = out_dispatcher_factory.create(art_params.vm, new InMemoryFastaFetch(), art_params.args);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.num_seqs() == 0) {
                        break;
                    }
                    job_id += 1;
                    SimulationJob sj(new InMemoryFastaFetch(std::move(fa_view)), coverage_info, job_id, true);
                    ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
                    job_pool.add(std::move(aje));
                }
                fasta_stream.close();
            }
        } else if (art_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
            if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                std::ifstream input_file_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream);
                auto [fasta_fetch, coverage_info] = batcher.fetch();
                input_file_stream.close();
                out_dispatcher = out_dispatcher_factory.create(art_params.vm, &fasta_fetch, art_params.args);

                InMemoryFastaStreamBatcher fsb(
                    static_cast<int>(fasta_fetch.num_seqs() / art_params.parallel + 1), &fasta_fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.num_seqs() == 0) {
                        break;
                    }
                    job_id += 1;
                    SimulationJob sj(new InMemoryFastaFetch(std::move(fa_view)), coverage_info, job_id, true);
                    ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
                    job_pool.add(std::move(aje));
                }
            } else { // Stream
                out_dispatcher
                    = out_dispatcher_factory.create(art_params.vm, new InMemoryFastaFetch(), art_params.args);
                std::ifstream pbsim3_transcript_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_params.batch_size, pbsim3_transcript_stream);
                while (true) {
                    auto [fa_view, coverage_info] = fsb.fetch();
                    if (fa_view.num_seqs() == 0) {
                        break;
                    }
                    job_id += 1;
                    SimulationJob sj(new InMemoryFastaFetch(std::move(fa_view)), coverage_info, job_id, true);
                    ArtJobExecutor aje(std::move(sj), art_params, out_dispatcher);
                    job_pool.add(std::move(aje));
                }
                pbsim3_transcript_stream.close();
            }
        } else {
            throw std::runtime_error("Unsupported input file type");
        }
        job_pool.stop();
        BOOST_LOG_TRIVIAL(info) << "Job pool stopped";
        out_dispatcher->close();
        delete out_dispatcher;
    }
}
} // namespace labw::art_modern
