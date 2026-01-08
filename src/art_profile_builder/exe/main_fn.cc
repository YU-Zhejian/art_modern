/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "art_modern_config.h" // For USED_HTSLIB_NAME

#include "art_profile_builder/exe/main_fn.hh"

#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/APBConstants.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/si_utils.hh"

#ifdef WITH_NCBI_NGS
#include <ncbi-vdb/NGS.hpp>

#include <ngs/ReadCollection.hpp>
#include <ngs/itf/ErrorMsg.hpp>
#endif

#include <concurrentqueue.h>

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "libam_support/utils/arithmetic_utils.hh"

namespace labw::art_modern {
constexpr std::size_t REPORT_SIZE = 10000000; // 10 Million reads
constexpr std::size_t APB_BATCH_SIZE = K_SIZE;

class ViewSamMtChild {

public:
    ViewSamMtChild(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
        const std::size_t thread_id, APBConfig config, moodycamel::ConcurrentQueue<std::vector<bam1_t*>>& read_queue)
        : ied1_(ied1)
        , ied2_(ied2)
        , thread_id_(thread_id)
        , config_(std::move(config))
        , read_queue_(read_queue)
    {
    }
    DELETE_COPY(ViewSamMtChild)
    DELETE_MOVE(ViewSamMtChild)
    void start() { worker_thread_ = std::thread(&ViewSamMtChild::run_, this); }
    ~ViewSamMtChild()
    {
        stop_();
        join();
    }
    void join()
    {
        if (worker_thread_.joinable()) {
            worker_thread_.join();
        }
    }

private:
    std::shared_ptr<IntermediateEmpDist> ied1_;
    std::shared_ptr<IntermediateEmpDist> ied2_;
    std::size_t thread_id_;
    APBConfig config_;
    moodycamel::ConcurrentQueue<std::vector<bam1_t*>>& read_queue_;
    std::atomic<bool> should_stop_ { false };
    std::condition_variable should_stop_cv_;
    std::mutex should_stop_mutex_;
    std::thread worker_thread_;
    am_readnum_t num_valid_reads_ = 0;
    am_readnum_t num_total_reads_ = 0;
    am_readnum_t num_total_in_ = 0;
    am_readnum_t num_wait_in_ = 0;

    void run_()
    {
        while (!should_stop_) {
            std::vector<bam1_t*> reads;
            if (read_queue_.try_dequeue(reads)) {
                if (reads.empty()) {
                    stop_();
                    break;
                }
                num_total_in_++;
                for (auto* b : reads) {
                    if (((config_.is_pe && ((b->core.flag & BAM_FREAD2) != 0)) ? ied2_ : ied1_)->parse_read(b)) {
                        num_valid_reads_++;
                    }
                    num_total_reads_++;
                    if (num_total_reads_ % REPORT_SIZE == 0) {
                        log_();
                    }
                    bam_destroy1(b);
                }
            } else {
                num_wait_in_++;
                std::unique_lock lk(should_stop_mutex_);
                should_stop_cv_.wait_for(lk, std::chrono::milliseconds(100), [this]() { return should_stop_.load(); });
            }
        }
        log_();

        BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id_
                                << ": Processed read batches: N. (IRetried/IAll): " << format_with_commas(num_wait_in_)
                                << "/" << format_with_commas(num_total_in_) << " ("
                                << (100.0 * static_cast<double>(num_wait_in_) / static_cast<double>(num_total_in_))
                                << "%).";
    }

    void stop_()
    {
        should_stop_.store(true);
        should_stop_cv_.notify_all();
    }
    void log_()
    {
        BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id_ << ": Processed "
                                << to_si(num_total_reads_, 2, static_cast<decltype(num_total_reads_)>(1000))
                                << " reads, "
                                << to_si(num_valid_reads_, 2, static_cast<decltype(num_total_reads_)>(1000)) << " ("
                                << static_cast<double>(num_valid_reads_) / static_cast<double>(num_total_reads_) * 100.0
                                << ")% valid reads.";
    }
};

void view_sam_mt(const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied1s,
    const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied2s, const std::size_t n_threads,
    const APBConfig& config)
{
    int retv = 0;

    htsThreadPool tpool = { nullptr, 0 };
    tpool.pool = CExceptionsProxy::assert_not_null(
        hts_tpool_init(static_cast<int>(config.num_io_threads)), USED_HTSLIB_NAME, "Failed to init HTS thread pool.");
    std::string modestr = "r";
    if (config.format == APB_FORMAT::BAM) {
        modestr += "b";
    } else if (config.format == APB_FORMAT::CRAM) {
        modestr += "c";
    } else if (config.format == APB_FORMAT::FASTQ) {
        modestr += "f";
    }
    auto* in = CExceptionsProxy::assert_not_null(
        hts_open(config.input_file_path.c_str(), modestr.c_str()), USED_HTSLIB_NAME, "Failed to open HTS file.");
    auto* hdr = CExceptionsProxy::assert_not_null(sam_hdr_read(in), USED_HTSLIB_NAME, "Failed to read SAM header.");
    auto* b = CExceptionsProxy::assert_not_null(bam_init1(), USED_HTSLIB_NAME, "Failed to init BAM record.");

    moodycamel::ConcurrentQueue<std::vector<bam1_t*>> read_queue { config.queue_size, 1, 0 };
    moodycamel::ProducerToken const read_producer_token { read_queue };

    hts_set_opt(in, HTS_OPT_THREAD_POOL, &tpool);
    std::vector<std::unique_ptr<ViewSamMtChild>> workers;
    for (std::size_t thread_id = 0; thread_id < n_threads; ++thread_id) {
        workers.emplace_back(
            std::make_unique<ViewSamMtChild>(ied1s[thread_id], ied2s[thread_id], thread_id, config, read_queue));
        workers.back()->start();
    }
    am_readnum_t num_parsed_reads = 0;
    am_readnum_t num_pushed_blocks = 0;
    am_readnum_t num_push_waits = 0;

    // Prepare cache
    bool at_eof = false;
    while (!at_eof) {
        std::vector<bam1_t*> batch;
        batch.reserve(APB_BATCH_SIZE);
        for (std::size_t i = 0; i < APB_BATCH_SIZE; ++i) {
            batch.emplace_back(bam_init1());
            retv = sam_read1(in, hdr, batch.back());
            num_parsed_reads += 1;
            if (num_parsed_reads >= config.first_n_reads) {
                at_eof = true;
                break;
            }
            if (retv == -1 /** EOF **/) {
                BOOST_LOG_TRIVIAL(info) << "Reached EOF when parsing read " << num_parsed_reads << ".";
                bam_destroy1(batch.back());
                batch.pop_back();
                at_eof = true;
                break;
            }
            if (retv < -1) {
                BOOST_LOG_TRIVIAL(info) << "EXCEPT when parsing read " << num_parsed_reads << ".";
                sam_hdr_destroy(hdr);
                bam_destroy1(b);
                hts_close(in);
                hts_tpool_destroy(tpool.pool);
                abort_mpi();
            }
        }
        bool success = read_queue.try_enqueue(read_producer_token, std::move(batch));
        if (!success) {
            num_push_waits++;
        }
        while (!success) {
            success = read_queue.try_enqueue(read_producer_token, std::move(batch));
        }
        num_pushed_blocks++;
    }
    // Send stop signals
    for (std::size_t i = 0; i < n_threads; ++i) {
        std::vector<bam1_t*> batch {};
        bool success = read_queue.try_enqueue(read_producer_token, std::move(batch));
        if (!success) {
            num_push_waits++;
        }
        while (!success) {
            success = read_queue.try_enqueue(read_producer_token, std::move(batch));
        }
        num_pushed_blocks++;
    }
    // Wait for workers to finish
    for (auto& worker : workers) {
        worker->join();
    }
    BOOST_LOG_TRIVIAL(info) << "Total read blocks pushed: N. (IRetried/IAll): " << format_with_commas(num_push_waits)
                            << "/" << format_with_commas(num_pushed_blocks) << " ("
                            << (100.0 * static_cast<double>(num_push_waits) / static_cast<double>(num_pushed_blocks))
                            << "%).";
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_close(in);
    hts_tpool_destroy(tpool.pool);
}
#ifdef WITH_NCBI_NGS

class ViewSraMtChild {

public:
    ViewSraMtChild(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
        const std::size_t thread_id, APBConfig config, const std::uint64_t start, const std::uint64_t end)
        : ied1_(ied1)
        , ied2_(ied2)
        , thread_id_(thread_id)
        , config_(std::move(config))
        , start_(start)
        , end_(end)
    {
    }
    DELETE_COPY(ViewSraMtChild)
    DELETE_MOVE(ViewSraMtChild)
    void start() { worker_thread_ = std::thread(&ViewSraMtChild::run_, this); }

    ~ViewSraMtChild() { join(); }
    void join()
    {
        if (worker_thread_.joinable()) {
            worker_thread_.join();
        }
    }

private:
    std::shared_ptr<IntermediateEmpDist> ied1_;
    std::shared_ptr<IntermediateEmpDist> ied2_;
    std::size_t thread_id_;
    APBConfig config_;
    std::uint64_t start_;
    std::uint64_t end_;

    std::thread worker_thread_;
    am_readnum_t num_valid_reads_ = 0;
    am_readnum_t num_total_reads_ = 0;

    void run_()
    {
        try {
            auto rcngs = ncbi::NGS::openReadCollection(config_.input_file_path);
            auto reads = rcngs.getReadRange(start_ + 1, end_ - start_);
            while (reads.nextRead()) {
                auto nfrags = reads.getNumFragments();
                num_total_reads_ += nfrags;
                if (config_.is_pe) {
                    if (nfrags != 2) {
                        continue;
                    }
                    reads.nextFragment();
                    auto base1 = reads.getFragmentBases().toString();
                    auto qual1 = reads.getFragmentQualities().toString();
                    reads.nextFragment();
                    auto base2 = reads.getFragmentBases().toString();
                    auto qual2 = reads.getFragmentQualities().toString();
                    if (ied1_->parse_read_ngs_fragment(base1, qual1) && ied2_->parse_read_ngs_fragment(base2, qual2)) {
                        num_valid_reads_ += 2;
                    }
                } else {
                    if (nfrags != 1) {
                        continue;
                    }
                    auto base = reads.getReadBases().toString();
                    auto qual = reads.getReadQualities().toString();
                    if (ied1_->parse_read_ngs_fragment(base, qual)) {
                        num_valid_reads_++;
                    }
                }

                if (num_total_reads_ % REPORT_SIZE == 0) {
                    log_();
                }
            }
            log_();

        } catch (const ngs::ErrorMsg& e) {
            BOOST_LOG_TRIVIAL(error) << "NGS error in thread " << thread_id_ << ": " << e.what();
            abort_mpi();
        }
    }
    void log_()
    {
        BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id_ << ": Processed "
                                << to_si(num_total_reads_, 2, static_cast<decltype(num_total_reads_)>(1000))
                                << " reads, "
                                << to_si(num_valid_reads_, 2, static_cast<decltype(num_total_reads_)>(1000)) << " ("
                                << static_cast<double>(num_valid_reads_) / static_cast<double>(num_total_reads_) * 100.0
                                << ")% valid reads.";
    }
};

void view_sra_mt(const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied1s,
    const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied2s, const std::size_t n_threads,
    const APBConfig& config)
{
    try {
        auto rcngs = ncbi::NGS::openReadCollection(config.input_file_path);
        auto const nreads = rcngs.getReadCount();
        auto const rc = am_min(nreads, static_cast<std::uint64_t>(config.first_n_reads));
        // Split work
        std::uint64_t chunk_size = rc / n_threads;
        std::vector<std::unique_ptr<ViewSraMtChild>> workers;
        for (std::size_t thread_id = 0; thread_id < n_threads; ++thread_id) {
            workers.emplace_back(std::make_unique<ViewSraMtChild>(ied1s[thread_id], ied2s[thread_id], thread_id, config,
                thread_id * chunk_size, (thread_id == n_threads - 1) ? rc : (thread_id + 1) * chunk_size));
            workers.back()->start();
        }
        // Wait for workers to finish
        for (auto& worker : workers) {
            worker->join();
        }

    } catch (const ngs::ErrorMsg& e) {
        BOOST_LOG_TRIVIAL(error) << "NGS error when opening SRA file: " << e.what();
        abort_mpi();
    }
}
#endif

} // namespace labw::art_modern
