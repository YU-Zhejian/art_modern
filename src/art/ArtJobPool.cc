#include "ArtJobPool.hh"
#include <thread>

namespace labw::art_modern {

#if defined(USE_NOP_PARALLEL)

void ArtJobPool::stop() { tick_.stop(); }

void ArtJobPool::add(ArtJobExecutor aje)
{
    tick_.add(&this_aje);
    std::scoped_lock lock(mutex_);
    aje.execute();
}

ArtJobPool::ArtJobPool(const ArtParams&)
    : tick_(art_params)
{
    tick_.start();
}

#elif defined(USE_ASIO_PARALLEL)
ArtJobPool::ArtJobPool(const ArtParams& art_params)
    : pool_(art_params.parallel)
{
}

void ArtJobPool::add(ArtJobExecutor aje)
{
    boost::asio::post(pool_, [this_aje = std::move(aje)]() mutable {
        this_aje.execute();
    });
}

void ArtJobPool::stop()
{
    pool_.join();
}

#endif

}
