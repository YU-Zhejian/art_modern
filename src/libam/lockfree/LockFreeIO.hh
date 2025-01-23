#pragma once
#include "art_modern_config.h" // For USE_ASIO_PARALLEL

#include "libam/Constants.hh" // For USE_ASIO_PARALLEL
#include "libam/utils/class_macros_utils.hh"

#include <boost/log/trivial.hpp>

#if !defined(USE_NOP_PARALLEL)
#include <concurrentqueue.h>
#endif

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <mutex> // NOLINT
#include <thread>

namespace labw::art_modern {

#if !defined(USE_NOP_PARALLEL)
template <typename T> class LockFreeIO {
public:
    DELETE_MOVE(LockFreeIO)
    DELETE_COPY(LockFreeIO)

    static const int QUEUE_SIZE = M_SIZE;
    static const int BULK_SIZE = K_SIZE;
    LockFreeIO()
        : queue_(QUEUE_SIZE)
    {
    }

    virtual ~LockFreeIO() { stop(); };
    void push(T&& value)
    {
        while (!queue_.try_enqueue(std::move(value))) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        num_reads_in_++;
    }
    void start() { thread_ = std::thread(&LockFreeIO::run, this); }

    void stop()
    {
        should_stop_ = true;
        if (thread_.joinable()) {
            thread_.join();
        }
    }
    virtual void write(T value) = 0;

private:
    moodycamel::ConcurrentQueue<T> queue_;
    std::atomic<bool> should_stop_ = false;
    std::thread thread_;
    std::atomic<std::size_t> num_reads_in_ = 0;
    std::atomic<std::size_t> num_reads_out_ = 0;

    void run()
    {
        std::size_t pop_ret_cnt = 0;
        std::array<T, BULK_SIZE> retp_a;
        while (!should_stop_) {
            pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), BULK_SIZE);

            num_reads_out_ += pop_ret_cnt;
            while (pop_ret_cnt > 0) {
                pop_ret_cnt--;
                write(std::move(retp_a[pop_ret_cnt]));
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), BULK_SIZE);
        while (pop_ret_cnt > 0) {
            num_reads_out_ += pop_ret_cnt;
            while (pop_ret_cnt > 0) {
                pop_ret_cnt--;
                write(std::move(retp_a[pop_ret_cnt]));
            }
            pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), BULK_SIZE);
        }
        BOOST_LOG_TRIVIAL(info) << "LockFreeIO::run() finished, consuming " << num_reads_in_ << " reads and writes "
                                << num_reads_out_ << " reads";
    }
};

#else
template <typename T> class LockFreeIO {
public:
    LockFreeIO() = default;

    virtual ~LockFreeIO()
    {
        BOOST_LOG_TRIVIAL(info) << "LockFreeIO::run() finished, consuming " << num_reads_in_ << " reads and writes "
                                << num_reads_out_ << " reads";
    };
    void push(T&& value)
    {
        std::scoped_lock lock(mutex_);
        num_reads_in_++;
        write(std::move(value));
        num_reads_out_++;
    }
    void start()
    {
        // Do nothing!
    }

    void stop()
    {
        // Do nothing!
    }
    virtual void write(T value) = 0;

private:
    std::atomic<std::size_t> num_reads_in_ = 0;
    std::atomic<std::size_t> num_reads_out_ = 0;
    std::mutex mutex_;
};
#endif
} // namespace labw::art_modern
