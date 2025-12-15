/**
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

#include "art_modern_config.h" // NOLINT: For USE_ASIO_PARALLEL

#include "libam_support/jobs/JobExecutor.hh"
#include "libam_support/jobs/JobPool.hh"

#include "libam_support/jobs/Scheduler.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/post.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_BS_PARALLEL, USE_ASIO_PARALLEL"
#endif

#if !defined(USE_NOP_PARALLEL)
#include <chrono>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>
#endif

#include <cstddef>
#include <memory>

namespace labw::art_modern {

namespace {
    class JobPoolSupervisor : public Scheduler<std::chrono::milliseconds> {
    public:
        JobPoolSupervisor(JobPool& jp, const std::size_t supervisor_interval_ms)
            : Scheduler(
                  std::chrono::milliseconds(supervisor_interval_ms), std::chrono::milliseconds(supervisor_interval_ms))
            , jp_(jp)
        {
        }
        ~JobPoolSupervisor() override { stop(); }
        DELETE_MOVE(JobPoolSupervisor)
        DELETE_COPY(JobPoolSupervisor)
        void callback() override { jp_.prune_finished_jobs(); }

    private:
        JobPool& jp_;
    };
} // namespace

std::size_t JobPool::n_running_executors() const { return cached_n_running_executors_; }

void JobPool::prune_finished_jobs()
{
#if defined(USE_NOP_PARALLEL)
#else
    const std::scoped_lock op_lock(operation_mutex_);
    std::vector<std::shared_ptr<JobExecutor>> passed_ajes_;
    for (auto aje : executors_) {
        if (aje && aje->is_running()) {
            passed_ajes_.emplace_back(std::move(aje));
        }
    }
    executors_ = std::move(passed_ajes_);
    cached_n_running_executors_ = executors_.size();
    cached_n_running_executors_cv_.notify_all();
#endif
}

void JobPool::stop()
{
#if defined(USE_NOP_PARALLEL)
#else
#if defined(USE_BS_PARALLEL)
    pool_.wait();
#elif defined(USE_ASIO_PARALLEL)
    pool_.join();
#endif
    (static_cast<JobPoolSupervisor*>(supervisor_))->stop();
#endif
    prune_finished_jobs(); // Should clear all jobs.
#ifdef CEU_CM_IS_DEBUG
    if (cached_n_running_executors_ != 0) {
        BOOST_LOG_TRIVIAL(fatal) << "Number of cached executors are not 0 when exiting the job pool! "
        << "Actual: " << cached_n_running_executors_;
        abort_mpi();
    }
#endif
}

void JobPool::add(const std::shared_ptr<JobExecutor>& job_executor)
{
#if defined(USE_NOP_PARALLEL)
    cached_n_running_executors_ = 1;
    job_executor->operator()();
    cached_n_running_executors_ = 0;
#else
    // Lock for adding jobs
    std::unique_lock add_lock(add_mutex_);
    // Spin until there's a slot
    while (cached_n_running_executors_ >= pool_size_) {
        cached_n_running_executors_cv_.wait_for(
            add_lock,
            std::chrono::milliseconds(add_spin_waits_ms_),
            [this]() { return this->cached_n_running_executors_ < this->pool_size_; }
    );
    }
    // Start the job first
#if defined(USE_BS_PARALLEL)
    [[maybe_unused]] auto future = pool_.submit_task([this_aje = job_executor]() mutable { this_aje->operator()(); });
#elif defined(USE_ASIO_PARALLEL)
    post(pool_, [this_aje = job_executor]() mutable { this_aje->operator()(); });
#endif
    // Then add to the list
    {
        const std::scoped_lock op_lock(operation_mutex_);
        executors_.emplace_back(job_executor);
    }
#endif
}

JobPool::~JobPool()
{
    stop();
    // No idea how to free supervisor.
}

JobPool::JobPool(const std::size_t pool_size)
#if defined(USE_NOP_PARALLEL)
    {}
#else
    : pool_(pool_size)
    , pool_size_(pool_size)
{
    supervisor_ = new JobPoolSupervisor(*this, supervisor_interval_ms_);
    static_cast<JobPoolSupervisor*>(supervisor_)->start();
}
#endif

} // namespace labw::art_modern
