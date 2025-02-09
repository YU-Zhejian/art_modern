/**
 * Implementation proxy of some third-party job pools.
 */
#pragma once
#include "art_modern_config.h" // NOLINT: For USE_ASIO_PARALLEL

#include "libam_support/jobs/JobExecutor.hh"
#include "libam_support/utils/class_macros_utils.hh"

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/thread_pool.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_BS_PARALLEL, USE_ASIO_PARALLEL"
#endif

#include <cstddef>
#include <memory>
#include <mutex>
#include <vector>

namespace labw::art_modern {

/**
 * Job pool implementation.
 */
class JobPool {
public:
    DELETE_COPY(JobPool)
    DELETE_MOVE(JobPool)

    /**
     * Default constructor.
     * @param pool_size Number of jobs that can be execuetd in parallel.
     */
    explicit JobPool(std::size_t pool_size);

    /**
     * Destructor.
     */
    ~JobPool();

    /**
     * Add a job to the pool.
     *
     * @param aje Job to add.
     */
    void add(const std::shared_ptr<JobExecutor>& aje);

    /**
     * Stop the pool.
     */
    void stop();
    /**
     * Number of job instances that are still alive in the pool.
     * @return As described.
     */
    [[nodiscard]] std::size_t n_running_ajes() const;

private:
#if defined(USE_NOP_PARALLEL)
    // Do nothing!
#elif defined(USE_BS_PARALLEL)
    BS::thread_pool<BS::tp::none> pool_;
#elif defined(USE_ASIO_PARALLEL)
    boost::asio::thread_pool pool_;
#endif
    std::vector<std::shared_ptr<JobExecutor>> ajes_;
    std::size_t pool_size_ = 1;
    std::mutex mutex_;
};

} // namespace labw::art_modern
