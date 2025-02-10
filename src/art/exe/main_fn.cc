#include "art_modern_config.h"

#include "art/exe/main_fn.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/ArtIOParams.hh"
#include "art/lib/ArtJobExecutor.hh"
#include "art/lib/ArtParams.hh"

#include "libam_support/Constants.hh"
#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/jobs/JobPool.hh"
#include "libam_support/jobs/SimulationJob.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/out/OutputDispatcher.hh"
#include "libam_support/ref/batcher/FastaStreamBatcher.hh"
#include "libam_support/ref/batcher/InMemoryFastaBatcher.hh"
#include "libam_support/ref/batcher/Pbsim3TranscriptBatcher.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/ref/fetch/FaidxFetch.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <atomic>
#include <chrono>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>

namespace labw::art_modern {

class JobPoolReporter {
public:
    explicit JobPoolReporter(const JobPool& jp)
        : jp_(jp)
    {
    }
    void stop()
    {
        should_stop_ = true;
        thread_.join();
    }
    void start() { thread_ = std::thread(&JobPoolReporter::job_, this); }

private:
    void job_()
    {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        while (!should_stop_) {
            BOOST_LOG_TRIVIAL(info) << "JobPoolReporter: " << jp_.n_running_ajes() << " JobExecutors running";
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }

    const JobPool& jp_;
    std::atomic<bool> should_stop_ { false };
    std::thread thread_;
};

class Generator {
public:
    const OutputDispatcherFactory out_dispatcher_factory;

    explicit Generator(const ArtParams& art_params, const ArtIOParams& art_io_params)
        : art_params_(art_params)
        , art_io_params_(art_io_params)
        , job_pool_(art_io_params.parallel)
        , reporter_(job_pool_)
    {
        reporter_.start();
    }
    void init_dispatcher(const std::shared_ptr<BaseFastaFetch>& fetch)
    {
        OutParams params{art_io_params_.parallel, art_io_params_.vm,  art_io_params_.args, fetch};
        out_dispatcher_ = out_dispatcher_factory.create(params);
    }

    void add(const std::shared_ptr<BaseFastaFetch>& fetch, const std::shared_ptr<CoverageInfo>& coverage_info)
    {
        SimulationJob sj { fetch, coverage_info, ++job_id_ };
        auto aje = std::make_shared<ArtJobExecutor>(std::move(sj), art_params_, out_dispatcher_);
        job_pool_.add(aje);
    }

    void wait()
    {
        BOOST_LOG_TRIVIAL(info) << "All jobs submitted. Waiting for job pool to stop...";
        job_pool_.stop();
        reporter_.stop();
        BOOST_LOG_TRIVIAL(info) << "Job pool stopped.";
        out_dispatcher_->close();
        BOOST_LOG_TRIVIAL(info) << "Output dispatchers cleared.";
    }

private:
    int job_id_ = 0;
    ArtParams art_params_;
    ArtIOParams art_io_params_;
    JobPool job_pool_;
    JobPoolReporter reporter_;
    std::shared_ptr<OutputDispatcher> out_dispatcher_;
};

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

void generate_all(const ArtParams& art_params, const ArtIOParams& art_io_params)
{
    Generator generator(art_params, art_io_params);
    if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
        const auto coverage_info
            = std::make_shared<CoverageInfo>(art_io_params.coverage_info.div(art_io_params.parallel));
        if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            std::shared_ptr<BaseFastaFetch> const fetch
                = std::make_shared<InMemoryFastaFetch>(art_io_params.input_file_name);
            generator.init_dispatcher(fetch);
            for (int i = 0; i < art_io_params.parallel; ++i) {
                generator.add(fetch, coverage_info);
            }
        } else {
            generator.init_dispatcher(std::make_shared<FaidxFetch>(art_io_params.input_file_name));
            for (int i = 0; i < art_io_params.parallel; ++i) {
                generator.add(std::make_shared<FaidxFetch>(art_io_params.input_file_name), coverage_info);
            }
        }
    } else {
        // Batch-based parallelism
        if (art_io_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {

            const auto coverage_info = std::make_shared<CoverageInfo>(art_io_params.coverage_info);
            if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                auto fetch = std::make_shared<InMemoryFastaFetch>(art_io_params.input_file_name);
                generator.init_dispatcher(fetch);

                InMemoryFastaBatcher fsb(static_cast<int>(fetch->num_seqs() / art_io_params.parallel + 1), fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
            } else {
                std::ifstream fasta_stream(art_io_params.input_file_name);
                FastaStreamBatcher fsb(art_io_params.batch_size, fasta_stream);
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
        } else if (art_io_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
            if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                std::ifstream input_file_stream(art_io_params.input_file_name);
                Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream);
                auto [fasta_fetch, coverage_info] = batcher.fetch();
                input_file_stream.close();
                generator.init_dispatcher(fasta_fetch);

                InMemoryFastaBatcher fsb(
                    static_cast<int>(fasta_fetch->num_seqs() / art_io_params.parallel + 1), fasta_fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
            } else { // Stream
                generator.init_dispatcher(std::make_shared<InMemoryFastaFetch>());
                std::ifstream pbsim3_transcript_stream(art_io_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_io_params.batch_size, pbsim3_transcript_stream);
                while (true) {
                    auto [fa_view, coverage_info] = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
                pbsim3_transcript_stream.close();
            }
        } else {
            throw std::runtime_error("Unsupported input file type");
        }
    }

    generator.wait();
    BOOST_LOG_TRIVIAL(info) << "Generator finished.";
}
} // namespace labw::art_modern
