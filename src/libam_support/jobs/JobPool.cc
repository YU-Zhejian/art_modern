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
#include "libam_support/utils/class_macros_utils.hh"

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/post.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_BS_PARALLEL, USE_ASIO_PARALLEL"
#endif

#include <atomic>
#include <chrono>
#include <cstddef>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace labw::art_modern {


std::size_t JobPool::n_running_executors() const
{
    return n_running_executors_;
}

void JobPool::prune_finished_jobs(){
    const std::scoped_lock op_lock(operation_mutex_);
    std::size_t n_running = 0;
    std::vector<std::shared_ptr<JobExecutor>> passed_ajes_;
    for (auto aje : executors_) {
        if (aje && aje->is_running()) {
            n_running++;
            passed_ajes_.emplace_back(std::move(aje));
        }
    }
    executors_ = std::move(passed_ajes_);
    n_running_executors_ = n_running;
}

void JobPool::stop()
{
#if defined(USE_NOP_PARALLEL)
    return;
#elif defined(USE_BS_PARALLEL)
    pool_.wait();
#elif defined(USE_ASIO_PARALLEL)
    pool_.join();
#endif
    should_stop_ = true;
    if (supervisor_thread_.joinable()) {
        supervisor_thread_.join();
    }
}

void JobPool::supervisor_() {
    while (! should_stop_){
        prune_finished_jobs();
        std::this_thread::sleep_for(std::chrono::milliseconds(supervisor_interval_ms_));
    }
}

void JobPool::add(const std::shared_ptr<JobExecutor>& job_executor)
{
    const std::scoped_lock add_lock(add_mutex_);
    // Spin until there's a slot
    while (n_running_executors_ >= pool_size_){
        prune_finished_jobs();
        std::this_thread::sleep_for(std::chrono::milliseconds(add_spin_waits_ms_));
    }
    // Start the job first
#if defined(USE_NOP_PARALLEL)
        aje->operator()();
#elif defined(USE_BS_PARALLEL)
        [[maybe_unused]] auto future = pool_.submit_task([this_aje = aje]() mutable { this_aje->operator()(); });
#elif defined(USE_ASIO_PARALLEL)
        post(pool_, [this_aje = job_executor]() mutable { this_aje->operator()(); });
#endif
    // Then add to the list
    {
        const std::scoped_lock op_lock(operation_mutex_);
        executors_.emplace_back(job_executor);
    }
}

JobPool::~JobPool() { stop(); }

JobPool::JobPool([[maybe_unused]] const std::size_t pool_size)
#if defined(USE_NOP_PARALLEL)
    {
    // No need to start the supervisor thread
    }
#elif defined(USE_BS_PARALLEL)
    : pool_(pool_size)
    , pool_size_(pool_size) { }
#elif defined(USE_ASIO_PARALLEL)
    : pool_(pool_size)
    , pool_size_(pool_size)
{
    supervisor_thread_ = std::thread(&JobPool::supervisor_, this);
}
#endif

} // namespace labw::art_modern
