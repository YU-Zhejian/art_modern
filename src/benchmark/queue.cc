/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "benchmark_utils.hh"

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

using namespace labw::art_modern; // NOLINT
int main()
{
    constexpr std::size_t NUM_TRIALS = 20;
    std::vector<std::size_t> times;
    for (std::size_t i = 0; i < NUM_TRIALS; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        bench_moody_camel_implicit();
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    }
    BOOST_LOG_TRIVIAL(info) << "Moody Camel (Implicit): " << describe(times) << " ms";
    times.clear();
    for (std::size_t i = 0; i < NUM_TRIALS; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        bench_moody_camel_explicit();
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    }
    BOOST_LOG_TRIVIAL(info) << "Moody Camel (Explicit): " << describe(times) << " ms";

    return EXIT_SUCCESS;
}
