#pragma once
#include "ArtJobExecutor.hh"

#if defined(PARALLEL_DISABLED)
#include <mutex>
#elif defined(PARALLEL_ASIO)
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#else
#error "No parallel strategy defined! One of: PARALLEL_DISABLED, PARALLEL_ASIO"
#endif

namespace labw {
namespace art_modern {
    class ArtJobPool {
    public:
        explicit ArtJobPool(const ArtParams& art_params);
        void add(ArtJobExecutor& aje);
        void stop();

    private:
#if defined(PARALLEL_DISABLED)
        std::mutex mutex_;
#elif defined(PARALLEL_ASIO)
        boost::asio::thread_pool pool_;
#endif
    };

}
}