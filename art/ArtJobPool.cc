#include "ArtJobPool.hh"

namespace labw::art_modern {
#if defined(PARALLEL_DISABLED)

void ArtJobPool::stop()
{
    // Do nothing!
}

void ArtJobPool::add(ArtJobExecutor aje)
{
    std::scoped_lock lock(mutex_);
    aje.execute();
}

ArtJobPool::ArtJobPool(const ArtParams&)
{
    // Do nothing!
}

#elif defined(PARALLEL_ASIO)
ArtJobPool::ArtJobPool(const ArtParams& art_params)
    : pool_(art_params.parallel)
{
}

void ArtJobPool::add(ArtJobExecutor aje)
{
    boost::asio::post(pool_, [this_aje = std::move(aje)]() mutable { this_aje.execute(); });
}

void ArtJobPool::stop() { pool_.join(); }

#endif

}
