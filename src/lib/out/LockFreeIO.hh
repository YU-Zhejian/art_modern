#pragma once

#include <boost/lockfree/queue.hpp>
#include <thread>
#include <random>

namespace labw::art_modern {

template <typename T> class LockFreeIO {
public:
    static const int QUEUE_SIZE = 65534;
    LockFreeIO()
        : queue(), dis_(100, 1000) {};
    ~LockFreeIO() = default;
    void push(T* value)
    {
        while (!queue.push(value)) {
            std::this_thread::sleep_for(std::chrono::nanoseconds (dis_(gen_)));
        }
    }
    void start() { thread_ = std::thread(&LockFreeIO::run, this); }

    void stop()
    {
        should_stop_ = true;
        if (thread_.joinable()) {
            thread_.join();
        }
    }
    virtual void write(T*) = 0;

private:
    boost::lockfree::queue<T*, boost::lockfree::fixed_sized<true>, boost::lockfree::capacity<QUEUE_SIZE>> queue;
    std::thread thread_;
    std::atomic<bool> should_stop_ = false;
    std::minstd_rand gen_;
    std::uniform_int_distribution<int> dis_;

    void run()
    {
        bool pop_ret = false;
        T* retp;
        while (!should_stop_) {
            pop_ret = queue.pop(retp);
            if (pop_ret) {
                write(retp);
            } else {
                // std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        pop_ret = queue.pop(retp);
        while (pop_ret) {
            write(retp);
            pop_ret = queue.pop(retp);
        }
    }
};

}