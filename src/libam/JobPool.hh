/**
 * Implementation proxy of some third-party job pools.
 */
#pragma once
#include "art_modern_config.h" // NOLINT: For USE_ASIO_PARALLEL

#include "libam/utils/class_macros_utils.hh"

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_BS_PARALLEL, USE_ASIO_PARALLEL"
#endif

#include <chrono>
#include <cstddef>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace labw::art_modern {

/**
 * Job pool implementation.
 * @tparam T Job that implements `operator()()`, with `is_running` property,
 * and can be held by `std::shared_ptr`.
 */
template <typename T> class JobPool {
public:
    DELETE_COPY(JobPool)
    DELETE_MOVE(JobPool)

    /**
     * Default constructor.
     * @param pool_size Number of jobs that can be execuetd in parallel.
     */
    explicit JobPool(const std::size_t pool_size)
#if defined(USE_NOP_PARALLEL)
        {}
#elif defined(USE_BS_PARALLEL)
        : pool_(pool_size)
        , pool_size_(pool_size) {}
#elif defined(USE_ASIO_PARALLEL)
        : pool_(pool_size)
        , pool_size_(pool_size)
    {
    }
#endif

    /**
     * Destructor.
     */
    ~JobPool()
    {
        stop();
    }

    /**
     * Add a job to the pool.
     *
     * @param aje Job to add.
     */
    void add(const std::shared_ptr<T>& aje)
    {
        const std::scoped_lock lock(mutex_);
        // Spin until there's a slot
        while (n_running_ajes() >= pool_size_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        ajes_.emplace_back(aje);
#if defined(USE_NOP_PARALLEL)
        aje->operator()();
#elif defined(USE_BS_PARALLEL)
        pool_.submit_task([this_aje = aje]() mutable { this_aje->operator()(); });
#elif defined(USE_ASIO_PARALLEL)
        post(pool_, [this_aje = aje]() mutable { this_aje->operator()(); });
#endif
    }

    /**
     * Stop the pool.
     */
    void stop()
    {
#if defined(USE_NOP_PARALLEL)
        return;
#elif defined(USE_BS_PARALLEL)
        pool_.wait();
#elif defined(USE_ASIO_PARALLEL)
        pool_.join();
#endif
    }
    /**
     * Number of job instances that are still alive in the pool.
     * @return As described.
     */
    [[nodiscard]] std::size_t n_running_ajes() const
    {
        std::size_t n_running = 0;
        for (const auto& aje : ajes_) {
            if (aje && aje->is_running) {
                n_running++;
            }
        }
        return n_running;
    }

private:
#if defined(USE_NOP_PARALLEL)
    std::mutex mutex_;
#elif defined(USE_BS_PARALLEL)
    BS::thread_pool<BS::tp::none> pool_;
    std::mutex mutex_;
#elif defined(USE_ASIO_PARALLEL)
    boost::asio::thread_pool pool_;
    std::mutex mutex_;
#endif
    std::vector<std::shared_ptr<T>> ajes_;
    std::size_t pool_size_ = 1;
};

} // namespace labw::art_modern
