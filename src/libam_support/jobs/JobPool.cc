#include "art_modern_config.h" // NOLINT: For USE_ASIO_PARALLEL

#include "libam_support/jobs/JobExecutor.hh"
#include "libam_support/jobs/JobPool.hh"

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/post.hpp>
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

std::size_t JobPool::n_running_ajes() const
{
    std::size_t n_running = 0;
    for (const auto& aje : ajes_) {
        if (aje && aje->is_running()) {
            n_running++;
        }
    }
    return n_running;
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
}

void JobPool::add(const std::shared_ptr<JobExecutor>& aje)
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
    [[maybe_unused]] auto future = pool_.submit_task([this_aje = aje]() mutable { this_aje->operator()(); });
#elif defined(USE_ASIO_PARALLEL)
    post(pool_, [this_aje = aje]() mutable { this_aje->operator()(); });
#endif
}

JobPool::~JobPool() { stop(); }

JobPool::JobPool(std::size_t pool_size)
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

} // namespace labw::art_modern
