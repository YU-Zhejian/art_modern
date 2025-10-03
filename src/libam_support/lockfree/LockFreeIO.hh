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

#pragma once

#include "absl/base/attributes.h"
#include "libam_support/Constants.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#include <concurrentqueue.h>

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <string>
#include <thread>
#include <utility>

namespace labw::art_modern {

template <typename T> class LockFreeIO {
public:
    DELETE_MOVE(LockFreeIO)
    DELETE_COPY(LockFreeIO)

    /**
     * Size of the moodycamel queue.
     */
    static constexpr int QUEUE_SIZE = M_SIZE;
    /**
     * Bulk size for dequeueing. Used for acceleration.
     */
    static constexpr int BULK_SIZE = K_SIZE;

    /** Sleep time for waiting if there's nothing to dequeue. **/
    static constexpr auto SLEEP_TIME = std::chrono::microseconds(10);

    /**
     * Constructor.
     * @param name Name of this IO worker, used for logging.
     */
    explicit LockFreeIO(std::string name);

    virtual ~LockFreeIO() = default;

    /**
     * Initialize the queue with the given number of producers.
     * @param num_explicit_producers Number of explicit producers (those with tokens).
     * @param num_implicit_producers Number of implicit producers (those without).
     */
    void init_queue(std::size_t num_explicit_producers, std::size_t num_implicit_producers);

    /**
     * Pushing a value into the queue.
     * @param value As described.
     */
    void push(T&& value);

    /**
     * Pushing a value into the queue with a producer token.
     * @param value As described.
     * @param token The producer token.
     */
    void push(T&& value, const ProducerToken& token);

    /**
     * Start the IO worker thread.
     */
    void start();
    /**
     * Stop the IO worker thread.
     */
    void stop();
    /**
     * Get a new producer token for this queue.
     * @return The producer token.
     */
    ProducerToken get_producer_token();

protected:
    /** Number of bytes written out. Should be set by each @link write @endlink  call **/
    std::size_t num_bytes_out_ = 0;
    /** Name of this IO worker, used for logging. **/
    const std::string name_;
    /**
     * Implementation of flushing and closing the actual filesystem, web socket, etc.
     * Default implementation does nothing.
     */
    virtual void flush_and_close();

    /**
     * Implementation of writing a value to the actual filesystem, web socket, etc.
     */
    virtual void write(T /**value**/);

private:
    /** Time point when the IO worker starts. **/
    std::chrono::high_resolution_clock::time_point start_time_;
    /** Time point when the IO worker ends. **/
    std::chrono::high_resolution_clock::time_point end_time_;
    /**
     * Number of @link push @endlink calls.
     */
    std::atomic<std::size_t> num_reads_in_ = 0;
    /**
     * Number of retries in @link push @endlink calls, indicating queue full.
     */
    std::atomic<std::size_t> num_wait_in_ = 0;

    /**
     * Number of @link write @endlink calls.
     */
    std::size_t num_reads_out_ = 0;
    /**
     * Number of bulk dequeue calls that return less than @link BULK_SIZE @endlink, indicating queue not full.
     */
    std::size_t num_out_not_full_ = 0;
    /**
     * Number of bulk dequeue calls that return @link BULK_SIZE @endlink.
     */
    std::size_t num_out_full_ = 0;
    /**
     * Number of bulk dequeue calls that return nothing, indicating queue not full.
     */
    std::size_t num_out_empty_ = 0;

    /**
     * Whether the @link log_()@endlink function has been called.
     */
    std::atomic<bool> had_logged_ = false;

    /** The lock-free queue. **/
    moodycamel::ConcurrentQueue<T> queue_;
    /** Whether the IO worker should stop. **/
    std::atomic<bool> should_stop_ = false;
    /** The IO worker thread. **/
    std::thread thread_;

    /** The main loop of the IO worker thread. **/
    void run_();
    /** Log the statistics. **/
    void log_() const;
};

template <typename T>
LockFreeIO<T>::LockFreeIO(std::string name)
    : name_(std::move(name))
    , queue_(QUEUE_SIZE)
{
}

template <typename T>
void LockFreeIO<T>::init_queue(const std::size_t num_explicit_producers, const std::size_t num_implicit_producers)
{
    queue_ = moodycamel::ConcurrentQueue<T>(QUEUE_SIZE, num_explicit_producers, num_implicit_producers);
}

template <typename T> void LockFreeIO<T>::push(T&& value)
{
    bool success = queue_.try_enqueue(std::move(value));
    if (!success) {
        ++num_wait_in_;
    }
    while (!success) {
        std::this_thread::sleep_for(SLEEP_TIME);
        success = queue_.try_enqueue(std::move(value));
    }
    ++num_reads_in_;
}

template <typename T>
ABSL_ATTRIBUTE_ALWAYS_INLINE inline void LockFreeIO<T>::push(T&& value, const ProducerToken& token)
{
    bool success = queue_.try_enqueue(token.token, std::move(value));
    if (!success) {
        ++num_wait_in_;
    }
    while (!success) {
        std::this_thread::sleep_for(SLEEP_TIME);
        success = queue_.try_enqueue(token.token, std::move(value));
    }
    ++num_reads_in_;
}

template <typename T> void LockFreeIO<T>::start()
{
    start_time_ = std::chrono::high_resolution_clock::now();
    thread_ = std::thread(&LockFreeIO::run_, this);
}

template <typename T> void LockFreeIO<T>::stop()
{
    should_stop_ = true;
    if (thread_.joinable()) {
        thread_.join();
    }
    flush_and_close();
    end_time_ = std::chrono::high_resolution_clock::now();
    if (!had_logged_) {
        log_();
    }
    had_logged_ = true;
}

template <typename T> void LockFreeIO<T>::flush_and_close() { }

template <typename T> void LockFreeIO<T>::write(T /**value**/) { }

template <typename T> ProducerToken LockFreeIO<T>::get_producer_token()
{
    return ProducerToken { moodycamel::ProducerToken { queue_ } };
}

template <typename T> void LockFreeIO<T>::run_()
{
    moodycamel::ConsumerToken token(queue_);
    std::size_t pop_ret_cnt = 0;
    std::array<T, BULK_SIZE> retp_a;
    while (!should_stop_) {
        pop_ret_cnt = queue_.try_dequeue_bulk(token, retp_a.data(), BULK_SIZE);
        if (pop_ret_cnt == 0) {
            num_out_empty_++;
            std::this_thread::sleep_for(SLEEP_TIME);
            continue;
        }
        if (pop_ret_cnt < BULK_SIZE) {
            num_out_not_full_++;
        } else {
            num_out_full_++;
        }

        num_reads_out_ += pop_ret_cnt;
        while (pop_ret_cnt > 0) {
            pop_ret_cnt--;
            write(std::move(retp_a[pop_ret_cnt]));
        }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    pop_ret_cnt = queue_.try_dequeue_bulk(token, retp_a.data(), BULK_SIZE);
    while (pop_ret_cnt > 0) {
        num_reads_out_ += pop_ret_cnt;
        while (pop_ret_cnt > 0) {
            pop_ret_cnt--;
            write(std::move(retp_a[pop_ret_cnt]));
        }
        pop_ret_cnt = queue_.try_dequeue_bulk(token, retp_a.data(), BULK_SIZE);
    }
}

template <typename T> void LockFreeIO<T>::log_() const
{
    const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_).count();
    BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Finished, consuming " << format_with_commas(num_reads_in_)
                            << " reads and writes " << format_with_commas(num_reads_out_) << " reads.";
    const auto num_writes = num_out_full_ + num_out_not_full_ + num_out_empty_;
    if (num_reads_in_ == 0 && num_reads_out_ == 0) {
        // Avoid division by zero
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Nothing written.";
        return;
    }
    BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: N. (IRetried/IAll): " << format_with_commas(num_wait_in_) << "("
                            << (100.0 * num_wait_in_ / num_reads_in_) << "%).";
    BOOST_LOG_TRIVIAL(info) << name_
                            << " LockFreeIO: N. (OFull/ONotFull/OEmpty): "
                            // OFull
                            << format_with_commas(num_out_full_) << "(" << (100.0 * num_out_full_ / num_writes)
                            << "%) / "
                            // ONotFull
                            << format_with_commas(num_out_not_full_) << "(" << (100.0 * num_out_not_full_ / num_writes)
                            << "%) / "
                            // OEmpty
                            << format_with_commas(num_out_empty_) << "(" << (100.0 * num_out_empty_ / num_writes)
                            << "%); "
                            // Avg. write batch size
                            << " avg. write batch size: " << (1.0 * num_reads_out_ / num_writes) << ".";
    BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: " << to_si(num_bytes_out_) << "B written in " << time / 1000.0
                            << " seconds. Speed: " << to_si(1.0 * num_bytes_out_ / (time / 1000.0)) << "B/s.";
}
} // namespace labw::art_modern
