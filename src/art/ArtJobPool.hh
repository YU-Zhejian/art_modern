#pragma once
#include "art_modern_config.h" // NOLINT: For USE_ASIO_PARALLEL

#include "art/ArtJobExecutor.hh"
#include "art/ArtParams.hh"

#if defined(USE_NOP_PARALLEL)
// Do nothing!
#elif defined(USE_BS_PARALLEL)
#include <BS_thread_pool.hpp>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio/thread_pool.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_ASIO_PARALLEL"
#endif

#include <cstddef>
#include <memory>
#include <mutex>
#include <vector>

namespace labw::art_modern {

class ArtJobPool {
public:
    explicit ArtJobPool(const ArtParams& art_params);
    void add(const std::shared_ptr<ArtJobExecutor>& aje);
    void stop();
    std::size_t n_running_ajes();

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
    std::vector<std::shared_ptr<ArtJobExecutor>> ajes_;
    std::size_t pool_size_ = 1;
};

} // namespace labw::art_modern
