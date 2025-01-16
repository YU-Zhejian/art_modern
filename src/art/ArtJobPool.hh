#pragma once
#include "art/ArtJobExecutor.hh"

#if defined(USE_NOP_PARALLEL)
#include <mutex>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio.hpp>
#include <mutex>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_ASIO_PARALLEL"
#endif

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
#elif defined(USE_ASIO_PARALLEL)
    boost::asio::thread_pool pool_;
    std::mutex mutex_;
#endif
    std::vector<std::shared_ptr<ArtJobExecutor>> ajes_;
    std::size_t pool_size_ = 1;
};

} // namespace labw::art_modern
