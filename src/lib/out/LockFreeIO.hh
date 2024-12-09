#pragma once

#include <concurrentqueue.h>
#include <thread>
#include <boost/log/trivial.hpp>

namespace labw::art_modern {

template <typename T> class LockFreeIO {
public:
    static const int QUEUE_SIZE = 1<<20;
    LockFreeIO(): queue_(QUEUE_SIZE){

    }

    virtual ~LockFreeIO() { stop(); };
    void push(std::unique_ptr<T>&& value)
    {
        while (!queue_.try_enqueue(std::move(value))) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        num_reads_in_ ++;
    }
    void start() { thread_ = std::thread(&LockFreeIO::run, this); }

    void stop()
    {
        should_stop_ = true;
        if (thread_.joinable()) {
            thread_.join();
        }
    }
    virtual void write(std::unique_ptr<T> value) = 0;

private:
    moodycamel::ConcurrentQueue<std::unique_ptr<T>> queue_;
    std::atomic<bool> should_stop_ = false;
    std::thread thread_;
    std::atomic<std::size_t > num_reads_in_ = 0;
    std::atomic<std::size_t > num_reads_out_ = 0;

    void run()
    {
        int pop_ret_cnt = 0;
        auto retp_a = std::array<std::unique_ptr<T>, 1<<10>();
        while (!should_stop_) {
            pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), 1<<10) - 1;
            while (pop_ret_cnt> 0) {
                write(std::move(retp_a[pop_ret_cnt]));
                pop_ret_cnt --;
                num_reads_out_++;
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), 1<<10) - 1;
        while (pop_ret_cnt> 0) {
            while (pop_ret_cnt> 0) {
                write(std::move(retp_a[pop_ret_cnt]));
                pop_ret_cnt --;
                num_reads_out_++;
            }
            pop_ret_cnt = queue_.try_dequeue_bulk(retp_a.data(), 1<<10) - 1;
        }
        BOOST_LOG_TRIVIAL(info) << "LockFreeIO::run() finished, consuming " << num_reads_in_ << " reads and writes " << num_reads_out_ << " reads";
    }
};

}