#include "art_modern_config.h"

#include "art/main_fn.hh"

#include "art/ArtConstants.hh"
#include "art/ArtJobExecutor.hh"
#include "art/ArtParams.hh"

#include "libam/Constants.hh"
#include "libam/JobPool.hh"
#include "libam/ds/CoverageInfo.hh"
#include "libam/jobs/SimulationJob.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/OutputDispatcher.hh"
#include "libam/ref/batcher/FastaStreamBatcher.hh"
#include "libam/ref/batcher/InMemoryFastaBatcher.hh"
#include "libam/ref/batcher/Pbsim3TranscriptBatcher.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/ref/fetch/FaidxFetch.hh"
#include "libam/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace labw::art_modern {

struct Generator {
    const OutputDispatcherFactory out_dispatcher_factory;

    explicit Generator(const ArtParams& art_params, const bool free_fasta_fetch_after_execution)
        : art_params_(art_params)
        , job_pool_(art_params.parallel)
        , free_fasta_fetch_after_execution_(free_fasta_fetch_after_execution)
    {
    }
    void init_dispatcher(const std::shared_ptr<BaseFastaFetch>& fetch)
    {
        out_dispatcher_ = out_dispatcher_factory.create(art_params_.vm, fetch.get(), art_params_.args);
    }

    void add(const std::shared_ptr<BaseFastaFetch>& fetch, const std::shared_ptr<CoverageInfo>& coverage_info)
    {
        SimulationJob sj { fetch, coverage_info, ++job_id_, free_fasta_fetch_after_execution_ };
        auto aje = std::make_shared<ArtJobExecutor>(std::move(sj), art_params_, out_dispatcher_);
        job_pool_.add(aje);
    }

    void wait()
    {
        job_pool_.stop();
        BOOST_LOG_TRIVIAL(info) << "Job pool stopped";
        out_dispatcher_->close();
    }

private:
    int job_id_ = 0;
    ArtParams art_params_;
    JobPool<ArtJobExecutor> job_pool_;
    const bool free_fasta_fetch_after_execution_;
    std::shared_ptr<BaseReadOutput> out_dispatcher_;
};

namespace {
    void generate_wgs(const ArtParams& art_params)
    {
        Generator generator(art_params, art_params.art_input_file_parser != INPUT_FILE_PARSER::MEMORY);

        const auto coverage_info = std::make_shared<CoverageInfo>(art_params.coverage_info.div(art_params.parallel));
        if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            std::shared_ptr<BaseFastaFetch> const fetch
                = std::make_shared<InMemoryFastaFetch>(art_params.input_file_name);
            generator.init_dispatcher(fetch);
            for (int i = 0; i < art_params.parallel; ++i) {
                generator.add(fetch, coverage_info);
            }
        } else {
            generator.init_dispatcher(std::make_shared<FaidxFetch>(art_params.input_file_name));
            for (int i = 0; i < art_params.parallel; ++i) {
                generator.add(std::make_shared<FaidxFetch>(art_params.input_file_name), coverage_info);
            }
        }
        generator.wait();
    }

} // namespace

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

void generate_all(const ArtParams& art_params)
{
    if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
        generate_wgs(art_params);
    } else {
        // Batch-based parallelism
        if (art_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {
            Generator generator(art_params, true);

            const auto coverage_info = std::make_shared<CoverageInfo>(art_params.coverage_info);
            if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                auto fetch = std::make_shared<InMemoryFastaFetch>(art_params.input_file_name);
                generator.init_dispatcher(fetch);

                InMemoryFastaBatcher fsb(static_cast<int>(fetch->num_seqs() / art_params.parallel + 1), fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
            } else {
                std::ifstream fasta_stream(art_params.input_file_name);
                FastaStreamBatcher fsb(art_params.batch_size, fasta_stream);
                generator.init_dispatcher(std::make_shared<InMemoryFastaFetch>());
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view.empty()) {
                        break;
                    }
                    generator.add(std::make_shared<InMemoryFastaFetch>(fa_view), coverage_info);
                }
                fasta_stream.close();
            }
            generator.wait();
        } else if (art_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
            Generator generator(art_params, true);
            if (art_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                std::ifstream input_file_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream);
                auto [fasta_fetch, coverage_info] = batcher.fetch();
                input_file_stream.close();
                generator.init_dispatcher(fasta_fetch);

                InMemoryFastaBatcher fsb(
                    static_cast<int>(fasta_fetch->num_seqs() / art_params.parallel + 1), fasta_fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
            } else { // Stream
                generator.init_dispatcher(std::make_shared<InMemoryFastaFetch>());
                std::ifstream pbsim3_transcript_stream(art_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_params.batch_size, pbsim3_transcript_stream);
                while (true) {
                    auto [fa_view, coverage_info] = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
                pbsim3_transcript_stream.close();
            }
            generator.wait();
        } else {
            throw std::runtime_error("Unsupported input file type");
        }
    }
}
} // namespace labw::art_modern
