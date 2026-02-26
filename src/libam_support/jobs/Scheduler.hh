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
#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace labw::art_modern {

template <typename Duration> class Scheduler {
public:
    Scheduler(Duration interval, Duration start_delay);

    virtual ~Scheduler();
    DELETE_COPY(Scheduler)
    DELETE_MOVE(Scheduler)

    virtual void callback() = 0;

    void start();
    void stop();

private:
    std::mutex mutex_;
    std::condition_variable condition_;
    std::atomic<bool> should_stop_ { false };
    Duration interval_;
    Duration start_delay_;
    std::thread thread_;
    void run_();
};

template <typename Duration> void Scheduler<Duration>::start() { thread_ = std::thread(&Scheduler::run_, this); }

template <typename Duration> void Scheduler<Duration>::stop()
{
    should_stop_.store(true);
    condition_.notify_all();
    if (thread_.joinable()) {
        thread_.join();
    }
}

template <typename Duration>
Scheduler<Duration>::Scheduler(const Duration interval, const Duration start_delay)
    : interval_(interval)
    , start_delay_(start_delay)
{
}

template <typename Duration> Scheduler<Duration>::~Scheduler() { stop(); }
template <typename Duration> void Scheduler<Duration>::run_()
{
    std::mutex mtx;
    std::unique_lock lock(mtx);
    condition_.wait_for(lock, start_delay_, [this]() { return this->should_stop_.load(); });

    while (!should_stop_.load()) {
        // Execute the callback
        callback();

        // Wait for 1 second or until stop signal
        condition_.wait_for(lock, interval_, [this]() { return this->should_stop_.load(); });
    }
}
} // namespace labw::art_modern
