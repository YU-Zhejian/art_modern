#include "libam_support/ds/PyQueue.hh"

#include "align_blkring.hh"

#include <boost/lockfree/queue.hpp> // NOLINT
#include <boost/log/trivial.hpp> // NOLINT

#include <concurrentqueue.h>

#include <array>
#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <string>
#include <thread>
#include <vector>

#define VERBOSE_IO 0 // NOLINT
#define FAST_RAND 1 // NOLINT
using namespace labw::art_modern; // NOLINT

namespace {

constexpr std::size_t data_len = 1U << 10U;
constexpr std::size_t queue_size = 1U << 14U;
constexpr std::size_t nitems = 1U << 20U;
constexpr std::size_t bulk_size = 1U << 5U;
constexpr std::size_t nthreads = 40;

char* randstr_cstr()
{
    auto* cstr = static_cast<char*>(std::malloc(sizeof(char) * data_len));
    return cstr;
}

std::string randstr()
{
    auto* cstr = randstr_cstr();
    std::string rets { cstr, data_len };
    std::free(cstr);
    return rets;
}

void pyqueue_producer(PyQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    {
        BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
    }
#endif
    std::size_t i = 0;
    while (i < nitems) {
        if (queue.put(randstr(), true)) {
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

void blkring_producer(RingBuffer<char*>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    {
        BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
    }
#endif
    std::size_t i = 0;
    while (i < nitems) {
        auto* string = randstr_cstr();
        if (queue.enqueue_ringbuf(&string) == 0) {
            i++;
        }
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " ended with " << i << " items enqueued";
#endif
}

void blkring_consumer(RingBuffer<char*>& queue)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread started";
#endif
    std::size_t i = 0;
    while (i < nitems * nthreads) {
        char** item = nullptr;
        if (queue.dequeue_ringbuf(&item) == 0) {
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

[[maybe_unused]] void mcqueue_explicit_producer(
    moodycamel::ConcurrentQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
#endif
    moodycamel::ProducerToken const token { queue };
    std::size_t i = 0;
    while (i < nitems) {
        while (!queue.try_enqueue(token, randstr())) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
#if VERBOSE_IO
            BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " halted with " << i
                                    << " items enqueued with qsize=" << queue.size_approx() << "/" << queue_size;
#endif
        }
        i++;
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " ended with " << i << " items enqueued";
#endif
}
[[maybe_unused]] void mcqueue_implicit_producer(
    moodycamel::ConcurrentQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " started";
#endif
    std::size_t i = 0;
    while (i < nitems) {
        while (!queue.try_enqueue(randstr())) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        i++;
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "producer thread " << id << " ended with " << i << " items enqueued";
#endif
}

[[maybe_unused]] void mcqueue_consumer_single(moodycamel::ConcurrentQueue<std::string>& queue)
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

void mcqueue_consumer_bulk_implicit(moodycamel::ConcurrentQueue<std::string>& queue)
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

void mcqueue_consumer_bulk_explicit(moodycamel::ConcurrentQueue<std::string>& queue)
{
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread started";
#endif
    moodycamel::ConsumerToken token(queue);
    std::size_t i = 0;
    std::size_t pop_ret_cnt = 0;
    std::array<std::string, bulk_size> retp_a;
    while (i < nitems * nthreads) {
        pop_ret_cnt = queue.try_dequeue_bulk(token, retp_a.data(), bulk_size);
#if VERBOSE_IO
        if (pop_ret_cnt == 0) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            BOOST_LOG_TRIVIAL(info) << i << "/" << nitems * nthreads << " elements consumed";
        }
#endif
        i += pop_ret_cnt;
    }
#if VERBOSE_IO
    BOOST_LOG_TRIVIAL(info) << "consumer thread ended";
#endif
}

void bench_blkring()
{
    RingBuffer<char*> rb { queue_size };
    std::vector<std::thread> producers;
    producers.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; i++) {
        producers.emplace_back(blkring_producer, std::ref(rb), i);
    }
    std::thread consumer(blkring_consumer, std::ref(rb));
    for (auto& producer : producers) {
        producer.join();
    }
    consumer.join();
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

void bench_moody_camel_explicit()
{
    moodycamel::ConcurrentQueue<std::string> queue(queue_size, nthreads, 0);
    std::vector<std::thread> producers;
    producers.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; i++) {
        producers.emplace_back(mcqueue_explicit_producer, std::ref(queue), i);
    }
    std::thread consumer(mcqueue_consumer_bulk_explicit, std::ref(queue));
    for (auto& producer : producers) {
        producer.join();
    }
    consumer.join();
}

void bench_moody_camel_implicit()
{
    moodycamel::ConcurrentQueue<std::string> queue(queue_size, 0, nthreads);
    std::vector<std::thread> producers;
    producers.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; i++) {
        producers.emplace_back(mcqueue_implicit_producer, std::ref(queue), i);
    }
    std::thread consumer(mcqueue_consumer_bulk_implicit, std::ref(queue));
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
    bench_blkring();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "RingBuffer: "
                            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms";

    start = std::chrono::high_resolution_clock::now();
    bench_pyqueue();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "PyQueue: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                            << " ms";

    start = std::chrono::high_resolution_clock::now();
    bench_moody_camel_implicit();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Moody Camel (Implicit): "
                            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms";
    start = std::chrono::high_resolution_clock::now();
    bench_moody_camel_explicit();
    end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Moody Camel (Explicit): "
                            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms";
    return EXIT_SUCCESS;
}
