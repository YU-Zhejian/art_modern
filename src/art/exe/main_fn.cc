/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "art_modern_config.h" // NOLINT

#include "art/exe/main_fn.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/ArtIOParams.hh"
#include "art/lib/ArtJobExecutor.hh"
#include "art/lib/ArtParams.hh"

#include "libam_support/Constants.hh"
#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ds/SkipLoaderSettings.hh"
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
#include "libam_support/utils/exception_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <thread>
#include <utility>

namespace labw::art_modern {

class JobPoolReporter {
public:
    explicit JobPoolReporter(JobPool& jp, const std::size_t reporting_interval_seconds = 1)
        : jp_(jp)
        , reporting_interval_seconds_(reporting_interval_seconds)
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
        std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
        while (!should_stop_) {
            BOOST_LOG_TRIVIAL(info) << "JobPoolReporter: " << jp_.n_running_ajes() << " JobExecutors running";
            std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
        }
    }

    JobPool& jp_;
    std::atomic<bool> should_stop_ { false };
    std::thread thread_;
    const std::size_t reporting_interval_seconds_;
};

class Generator {
public:
    const OutputDispatcherFactory out_dispatcher_factory;

    Generator(ArtParams art_params, const ArtIOParams& art_io_params)
        : art_params_(std::move(art_params))
        , art_io_params_(art_io_params)
        , job_pool_(art_io_params.parallel)
        , reporter_(job_pool_)
    {
        reporter_.start();
    }

    void init_dispatcher(const std::shared_ptr<BaseFastaFetch>& fetch, const bool clear_after_use)
    {
        clear_after_use_ = clear_after_use;
        OutParams const params { art_io_params_.parallel, art_io_params_.vm, art_io_params_.args, fetch };
        out_dispatcher_ = out_dispatcher_factory.create(params);
    }

    void add(const std::shared_ptr<BaseFastaFetch>& fetch, const std::shared_ptr<CoverageInfo>& coverage_info)
    {
        SimulationJob sj { fetch, coverage_info, ++job_id_ };
        const auto aje = std::make_shared<ArtJobExecutor>(std::move(sj), art_params_, out_dispatcher_, clear_after_use_);
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
    bool clear_after_use_ = true;
    std::size_t job_id_ = 0;
    const ArtParams art_params_;
    const ArtIOParams art_io_params_;
    JobPool job_pool_;
    JobPoolReporter reporter_;
    std::shared_ptr<OutputDispatcher> out_dispatcher_;
};

void print_banner()
{
    if (!is_on_mpi_main_process_or_nompi()) {
        return;
    }
    BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina (" << ART_PROGRAM_NAME << ") v. " ART_MODERN_VERSION
                            << " at <" << ART_MODERN_URL << ">";
    BOOST_LOG_TRIVIAL(info) << "Based on ART_Illumina: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
    BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
    BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <yuzj25@seas.upenn.edu>";
#ifdef CEU_CM_IS_DEBUG
    BOOST_LOG_TRIVIAL(info) << "Debugging functions enabled.";
#endif
}

void generate_all(const ArtParams& art_params, const ArtIOParams& art_io_params)
{
    Generator generator(art_params, art_io_params);
    if (art_params.art_simulation_mode == SIMULATION_MODE::WGS) {
        const std::size_t div_by = art_io_params.parallel * (have_mpi() ? mpi_size() : 1);
        const auto coverage_info = std::make_shared<CoverageInfo>(art_io_params.coverage_info.div(div_by));
        if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
            std::shared_ptr<BaseFastaFetch> const fetch
                = std::make_shared<InMemoryFastaFetch>(art_io_params.input_file_name);
            generator.init_dispatcher(fetch, false);
            for (std::size_t i = 0; i < art_io_params.parallel; ++i) {
                generator.add(fetch, coverage_info);
            }
        } else {
            generator.init_dispatcher(std::make_shared<FaidxFetch>(art_io_params.input_file_name), true);
            for (std::size_t i = 0; i < art_io_params.parallel; ++i) {
                generator.add(std::make_shared<FaidxFetch>(art_io_params.input_file_name), coverage_info);
            }
        }
    } else {
        // Batch-based parallelism
        const auto sls = SkipLoaderSettings::from_mpi();
        if (art_io_params.art_input_file_type == INPUT_FILE_TYPE::FASTA) {
            const auto coverage_info = std::make_shared<CoverageInfo>(art_io_params.coverage_info);
            std::ifstream fasta_stream(art_io_params.input_file_name);
            if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                FastaStreamBatcher fsb_f(std::numeric_limits<std::size_t>::max(), fasta_stream, sls);
                auto fetch = std::make_shared<InMemoryFastaFetch>(fsb_f.fetch());
                generator.init_dispatcher(fetch, true);
                InMemoryFastaBatcher fsb(static_cast<int>(fetch->num_seqs() / art_io_params.parallel + 1), fetch);
                while (true) {
                    auto fa_view = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
            } else {
                FastaStreamBatcher fsb_f(art_io_params.batch_size, fasta_stream, sls);
                generator.init_dispatcher(std::make_shared<InMemoryFastaFetch>(), true);
                while (true) {
                    auto fa_view = fsb_f.fetch();
                    if (fa_view.empty()) {
                        break;
                    }
                    generator.add(std::make_shared<InMemoryFastaFetch>(std::move(fa_view)), coverage_info);
                }
            }
            fasta_stream.close();
        } else if (art_io_params.art_input_file_type == INPUT_FILE_TYPE::PBSIM3_TRANSCRIPTS) {
            if (art_io_params.art_input_file_parser == INPUT_FILE_PARSER::MEMORY) {
                std::ifstream input_file_stream(art_io_params.input_file_name);
                Pbsim3TranscriptBatcher batcher(std::numeric_limits<int>::max(), input_file_stream, sls);
                const auto [fasta_fetch, coverage_info] = batcher.fetch();
                input_file_stream.close();
                generator.init_dispatcher(fasta_fetch, true);

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
                generator.init_dispatcher(std::make_shared<InMemoryFastaFetch>(), true);
                std::ifstream pbsim3_transcript_stream(art_io_params.input_file_name);
                Pbsim3TranscriptBatcher fsb(art_io_params.batch_size, pbsim3_transcript_stream, sls);
                while (true) {
                    const auto [fa_view, coverage_info] = fsb.fetch();
                    if (fa_view->empty()) {
                        break;
                    }
                    generator.add(fa_view, coverage_info);
                }
                pbsim3_transcript_stream.close();
            }
        }
    }

    generator.wait();
    BOOST_LOG_TRIVIAL(info) << "Generator finished.";
}
} // namespace labw::art_modern
