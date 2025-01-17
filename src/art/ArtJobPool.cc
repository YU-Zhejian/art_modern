#include "art/ArtJobPool.hh"

#include "art/ArtJobExecutor.hh"
#include "art/ArtParams.hh"

#include <boost/asio/post.hpp>

#include <chrono>
#include <cstddef>
#include <memory>
#include <mutex>
#include <thread>

namespace labw::art_modern {

#if defined(USE_NOP_PARALLEL)

void ArtJobPool::stop() { }

void ArtJobPool::add(const std::shared_ptr<ArtJobExecutor>& aje)
{
    std::scoped_lock lock(mutex_);
    ajes_.emplace_back(aje);
    aje->operator()();
}

ArtJobPool::ArtJobPool(const ArtParams&) { }

#elif defined(USE_ASIO_PARALLEL)
ArtJobPool::ArtJobPool(const ArtParams& art_params)
    : pool_(art_params.parallel)
    , pool_size_(art_params.parallel)
{
}

void ArtJobPool::add(const std::shared_ptr<ArtJobExecutor>& aje)
{
    std::scoped_lock lock(mutex_);
    // Spin until there's a slot
    while (n_running_ajes() >= pool_size_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    ajes_.emplace_back(aje);
    post(pool_, [this_aje = aje]() mutable { this_aje->operator()(); });
}

void ArtJobPool::stop() { pool_.join(); }

#endif

std::size_t ArtJobPool::n_running_ajes()
{
    std::size_t n_running = 0;
    for (const auto& aje : ajes_) {
        if (aje && aje->is_running) {
            n_running++;
        }
    }
    return n_running;
}
} // namespace labw::art_modern
