#include "libam/ds/PyQueue.hh"

#include <boost/lockfree/queue.hpp> // NOLINT
#include <boost/log/trivial.hpp> // NOLINT

#include <concurrentqueue.h>

#include <array>
#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>
#include <thread>
#include <vector>

#define VERBOSE_IO 0 // NOLINT
#define FAST_RAND 1 // NOLINT
using namespace labw::art_modern; // NOLINT

namespace {

constexpr std::size_t data_len = 1U << 10U;
constexpr std::size_t queue_size = 1U << 10U;
constexpr std::size_t nitems = 1U << 15U;
constexpr std::size_t bulk_size = 1U << 5U;
constexpr std::size_t nthreads = 40;

using randgen_t = std::minstd_rand;

std::string randstr([[maybe_unused]] randgen_t& gen)
{
#if FAST_RAND
    char* cstr = static_cast<char*>(std::malloc(sizeof(char) * data_len));
    std::string rets { cstr, data_len };
    std::free(cstr);
    return rets;
#else
    std::uniform_int_distribution<char> dist(CHAR_MIN, CHAR_MAX);
    std::string rets;
    rets.resize(strlen);
    for (auto& ch : rets) {
        ch = dist(gen);
    }
    return std::move(rets);
#endif
}

void pyqueue_producer(PyQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    {
        BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
    }
#endif
    randgen_t gen { std::random_device()() };
    std::size_t i = 0;
    while (i < nitems) {
        if (queue.put(randstr(gen), true)) {
            i++;
        }
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " ended with " << i << " items enqueued";
#endif
}

void pyqueue_consumer(PyQueue<std::string>& queue)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread started";
#endif
    std::string item;
    std::size_t i = 0;
    while (i < nitems * nthreads) {
        if (queue.get(item, false)) {
            i++;
        }
#if VERBOSE_IO
        else {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            BOOST_LOG_TRIVIAL(info) << i << " elements consumed";
        }
#endif
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread ended";
#endif
}

void mcqueue_producer(moodycamel::ConcurrentQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
#endif
    moodycamel::ProducerToken const token(queue);
    randgen_t gen { std::random_device()() };
    std::size_t i = 0;
    while (i < nitems) {
        while (!queue.try_enqueue(randstr(gen))) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        i++;
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " ended with " << i << " items enqueued";
#endif
}

void mcqueue_consumer_single(moodycamel::ConcurrentQueue<std::string>& queue)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread started";
#endif
    std::string item;
    std::size_t i = 0;
    bool pop_ret_cnt = false;
    while (i < nitems * nthreads) {
        pop_ret_cnt = queue.try_dequeue(item);
#if VERBOSE_IO
        if (!pop_ret_cnt) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            BOOST_LOG_TRIVIAL(info) << i << " elements consumed";
            BOOST_LOG_TRIVIAL(info) << "queue size: " << queue.size_approx();
        }
#endif
        if (pop_ret_cnt) {
            i++;
        }
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread ended";
#endif
}

void mcqueue_consumer_bulk(moodycamel::ConcurrentQueue<std::string>& queue)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread started";
#endif
    std::size_t i = 0;
    std::size_t pop_ret_cnt = 0;
    std::array<std::string, bulk_size> retp_a;
    while (i < nitems * nthreads) {
        pop_ret_cnt = queue.try_dequeue_bulk(retp_a.data(), bulk_size);
#if VERBOSE_IO
        if (pop_ret_cnt == 0) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            BOOST_LOG_TRIVIAL(info) << i << " elements consumed";
        }
#endif
        i += pop_ret_cnt;
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread ended";
#endif
}

void bench_pyqueue()
{
    PyQueue<std::string> queue(queue_size);
    std::vector<std::thread> producers;
    producers.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; i++) {
        producers.emplace_back(pyqueue_producer, std::ref(queue), i);
    }
    std::thread consumer(pyqueue_consumer, std::ref(queue));
    for (auto& producer : producers) {
        producer.join();
    }
    consumer.join();
}

void bench_moody_camel()
{
    moodycamel::ConcurrentQueue<std::string> queue(queue_size);
    std::vector<std::thread> producers;
    producers.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; i++) {
        producers.emplace_back(mcqueue_producer, std::ref(queue), i);
    }
    std::thread consumer(mcqueue_consumer_bulk, std::ref(queue));
    for (auto& producer : producers) {
        producer.join();
    }
    consumer.join();
}

} // namespace
int main()
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    start = std::chrono::high_resolution_clock::now();
    bench_pyqueue();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "PyQueue: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                            << " ms";
    start = std::chrono::high_resolution_clock::now();
    bench_moody_camel();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Moody Camel: "
                            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms";

    return EXIT_SUCCESS;
}
