#include "ArtJobPool.hh"

namespace labw {
namespace art_modern
{
#if defined(PARALLEL_DISABLED)

void ArtJobPool::stop()
{
// Do nothing!
}

void ArtJobPool::add(ArtJobExecutor aje)
{
    std::lock_guard<std::mutex> lock(mutex_);
    aje.execute();
}

ArtJobPool::ArtJobPool(const ArtParams &)
{
// Do nothing!
}

#elif defined(PARALLEL_ASIO)
ArtJobPool::ArtJobPool(const ArtParams &art_params)
:pool_(art_params.parallel == PARALLEL_ALL? static_cast<int>(std::thread::hardware_concurrency()) : art_params.parallel)
{

}

void ArtJobPool::add(ArtJobExecutor aje)
{
    auto func = [&aje] {
        aje.execute();
    };
    boost::asio::post(pool_, func);
}

void ArtJobPool::stop()
{
pool_.join();
}

#endif

}
}