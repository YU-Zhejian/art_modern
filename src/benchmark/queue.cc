#include "ds/PyQueue.hh"
#include <array>
#include <boost/lockfree/queue.hpp>
#include <chrono>
#include <climits>
#include <concurrentqueue.h>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

#define VERBOSE_IO 0

using namespace labw::art_modern;

constexpr std::size_t data_len = (1 << 10);
constexpr std::size_t queue_size = (1 << 10);
constexpr std::size_t nitems = 1 << 15;
constexpr std::size_t bulk_size = (1 << 5);
constexpr std::size_t nthreads = 40;

std::mutex iolock;

using randgen_t = std::minstd_rand;
#if (1)
std::string randstr([[maybe_unused]] randgen_t& gen)
{
    char* cstr = static_cast<char*>(std::malloc(sizeof(char) * data_len));
    std::string rets { cstr, data_len };
    std::free(cstr);
    return std::move(rets);
}
#else
std::string randstr(randgen_t& gen)
{
    std::uniform_int_distribution<char> dist(CHAR_MIN, CHAR_MAX);
    std::string rets;
    rets.resize(strlen);
    for (auto& ch : rets) {
        ch = dist(gen);
    }
    return std::move(rets);
}
#endif

void pyqueue_producer(PyQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "producer thread " << id << " started" << std::endl;
    }
#endif
    randgen_t gen;
    std::size_t i = 0;
    while (i < nitems) {
        std::string item = randstr(gen);
        if (queue.put(std::move(item), true)) {
            i++;
        }
    }
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "producer thread " << id << " ended with " << i << " items enqueued" << std::endl;
    }
#endif
}
void pyqueue_consumer(PyQueue<std::string>& queue)
{
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread started" << std::endl;
    }
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
            {
                std::scoped_lock lock(iolock);
                std::cout << i << " elements consumed" << std::endl;
            }
        }
#endif
    }
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread ended" << std::endl;
    }
#endif
}

void mcqueue_producer(moodycamel::ConcurrentQueue<std::string>& queue, [[maybe_unused]] std::size_t id)
{
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "producer thread " << id << " started" << std::endl;
    }
#endif
    moodycamel::ProducerToken token(queue);
    randgen_t gen;
    std::size_t i = 0;
    while (i < nitems) {
        std::string item = randstr(gen);
        while (!queue.try_enqueue(std::move(item))) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        i++;
    }
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "producer thread " << id << " ended with " << i << " items enqueued" << std::endl;
    }
#endif
}

void mcqueue_consumer_single(moodycamel::ConcurrentQueue<std::string>& queue)
{
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread started" << std::endl;
    }
#endif
    std::string item;
    std::size_t i = 0;
    bool pop_ret_cnt;
    while (i < nitems * nthreads) {
        pop_ret_cnt = queue.try_dequeue(item);
#if VERBOSE_IO
        if (!pop_ret_cnt) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            {
                std::scoped_lock lock(iolock);
                std::cout << i << " elements consumed" << std::endl;
                std::cout << "queue size: " << queue.size_approx() << std::endl;
            }
        }
#endif
        if (pop_ret_cnt) {
            i++;
        }
    }
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread ended" << std::endl;
    }
#endif
}

void mcqueue_consumer_bulk(moodycamel::ConcurrentQueue<std::string>& queue)
{
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread started" << std::endl;
    }
#endif
    std::string item;
    std::size_t i = 0;
    std::size_t pop_ret_cnt;
    std::array<std::string, bulk_size> retp_a;
    while (i < nitems * nthreads) {
        pop_ret_cnt = queue.try_dequeue_bulk(retp_a.data(), bulk_size);
#if VERBOSE_IO
        if (pop_ret_cnt == 0) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            {
                std::scoped_lock lock(iolock);
                std::cout << i << " elements consumed" << std::endl;
            }
        }
#endif
        i += pop_ret_cnt;
    }
#if VERBOSE_IO
    {
        std::scoped_lock lock(iolock);
        std::cout << "consumer thread ended" << std::endl;
    }
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

int main()
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    bench_pyqueue();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "PyQueue: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
    start = std::chrono::high_resolution_clock::now();
    bench_moody_camel();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Moody Camel: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;

    return EXIT_SUCCESS;
}
