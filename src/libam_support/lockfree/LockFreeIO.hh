#pragma once

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
#include <mutex> // NOLINT
#include <string>
#include <thread>
#include <utility>

namespace labw::art_modern {

template <typename T> class LockFreeIO {
public:
    DELETE_MOVE(LockFreeIO)
    DELETE_COPY(LockFreeIO)

    static const int QUEUE_SIZE = M_SIZE;
    static const int BULK_SIZE = K_SIZE;

    constexpr static std::chrono::duration sleep_time = std::chrono::microseconds(10);

    explicit LockFreeIO(std::string name)
        : name_(std::move(name))
        , queue_(QUEUE_SIZE)
    {
    }

    void init_queue(const std::size_t num_explicit_producers, const std::size_t num_implicit_producers)
    {
        queue_ = moodycamel::ConcurrentQueue<T>(QUEUE_SIZE, num_explicit_producers, num_implicit_producers);
    }

    virtual ~LockFreeIO() = default;

    void push(T&& value)
    {
        while (!queue_.try_enqueue(std::move(value))) {
            ++num_wait_in_;
            std::this_thread::sleep_for(sleep_time);
        }
        ++num_reads_in_;
    }

    void push(T&& value, const ProducerToken& token)
    {
        while (!queue_.try_enqueue(token.token, std::move(value))) {
            ++num_wait_in_;
            std::this_thread::sleep_for(sleep_time);
        }
        ++num_reads_in_;
    }

    void start()
    {
        start_time_ = std::chrono::high_resolution_clock::now();
        thread_ = std::thread(&LockFreeIO::run, this);
    }
    virtual void flush_and_close() {};

    void stop()
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

    virtual void write(T value) = 0;
    ProducerToken get_producer_token() { return ProducerToken { moodycamel::ProducerToken { queue_ } }; }

protected:
    std::size_t num_bytes_out_ = 0;
    const std::string name_;

private:
    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point end_time_;
    std::atomic<std::size_t> num_reads_in_ = 0;
    std::atomic<std::size_t> num_wait_in_ = 0;

    std::size_t num_reads_out_ = 0;
    std::size_t num_wait_out_not_full_ = 0;
    std::size_t num_wait_out_empty_ = 0;
    std::size_t num_nowait_out_ = 0;

    std::atomic<bool> had_logged_ = false;

    moodycamel::ConcurrentQueue<T> queue_;
    std::atomic<bool> should_stop_ = false;
    std::thread thread_;

    void run()
    {
        moodycamel::ConsumerToken token(queue_);
        std::size_t pop_ret_cnt = 0;
        std::array<T, BULK_SIZE> retp_a;
        while (!should_stop_) {
            pop_ret_cnt = queue_.try_dequeue_bulk(token, retp_a.data(), BULK_SIZE);
            if (pop_ret_cnt == 0) {
                num_wait_out_empty_++;
                std::this_thread::sleep_for(sleep_time);
                continue;
            }
            if (pop_ret_cnt < BULK_SIZE) {
                num_wait_out_not_full_++;
            } else {
                num_nowait_out_++;
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
    void log_() const
    {
        const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_).count();
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Finished, consuming " << format_with_commas(num_reads_in_)
                                << " reads and writes " << format_with_commas(num_reads_out_) << " reads.";
        BOOST_LOG_TRIVIAL(info)
            << name_ << " LockFreeIO: N. Waitings (I/ONotFull/OEmpty): " << format_with_commas(num_wait_in_) << " / "
            << format_with_commas(num_wait_out_not_full_) << "("
            << (100.0 * num_wait_out_not_full_ / (num_wait_out_empty_ + num_wait_out_not_full_ + num_nowait_out_))
            << "%) / " << format_with_commas(num_wait_out_empty_) << "("
            << (100.0 * num_wait_out_empty_ / (num_wait_out_empty_ + num_wait_out_not_full_ + num_nowait_out_))
            << "%).";
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: " << to_si(num_bytes_out_) << "B written in " << time / 1000.0
                                << " seconds. Speed: " << to_si(1.0 * num_bytes_out_ / (time / 1000.0)) << "B/s.";
    }
};

} // namespace labw::art_modern
