#pragma once
#include "ArtJobExecutor.hh"

#if defined(USE_NOP_PARALLEL)
#include <mutex>
#elif defined(USE_ASIO_PARALLEL)
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#else
#error "No parallel strategy defined! One of: USE_NOP_PARALLEL, USE_ASIO_PARALLEL"
#endif

namespace labw::art_modern {

    class ArtJobPool {
public:
    explicit ArtJobPool(const ArtParams& art_params);
    void add(ArtJobExecutor aje);
    void stop();

private:
#if defined(USE_NOP_PARALLEL)
    std::mutex mutex_;
#elif defined(USE_ASIO_PARALLEL)
    boost::asio::thread_pool pool_;
#endif
};

}
