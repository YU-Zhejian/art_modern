#include "ArtJobPool.hh"
#include <boost/log/trivial.hpp>
#include <thread>

namespace labw::art_modern {

Tick::Tick(const ArtParams& art_params)
    : read_size(art_params.read_len * (art_params.art_lib_const_mode == ART_LIB_CONST_MODE::SE ? 1 : 2))
{
}

void Tick::start() { thread_ = std::thread(&Tick::run, this); }

void Tick::stop()
{
    should_stop_ = true;
    if (thread_.joinable()) {
        thread_.join();
    }
}

void Tick::run()
{
    return; // FIXME: The following code have bugs. Will segfault on AJEs that are dead.
    std::chrono::time_point start = std::chrono::system_clock::now();
    size_t prev_num_reads = 0;
    int alive_ajes;
    size_t num_reads;
    while (!should_stop_) {
        std::this_thread::sleep_for(std::chrono::seconds(sleep_time));
        num_reads = 0;
        alive_ajes = 0;
        for (auto& aje : ajes_) {
            if (aje->is_running) {
                num_reads += aje->num_reads;
                alive_ajes++;
            }
        }
        std::chrono::time_point now = std::chrono::system_clock::now();
        const auto num_secs = std::chrono::duration_cast<std::chrono::seconds>((now - start)).count();
        BOOST_LOG_TRIVIAL(info) << "Tick: " << num_reads
                                << " reads generated. Speed: " << (num_reads - prev_num_reads) * read_size / sleep_time
                                << " nt/s; Avg. Speed: " << num_reads * read_size / num_secs << " nt/s; " << alive_ajes
                                << " AJEs alive";
        prev_num_reads = num_reads;
    }
}

void Tick::add(ArtJobExecutor* aje) { ajes_.emplace_back(aje); }

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
    , tick_(art_params)
{
    // tick_.start();
}

void ArtJobPool::add(ArtJobExecutor aje)
{
    boost::asio::post(pool_, [this_aje = std::move(aje), this]() mutable {
        // tick_.add(&this_aje);
        this_aje.execute();
    });
}

void ArtJobPool::stop()
{
    pool_.join();
    // tick_.stop();
}

#endif

}
