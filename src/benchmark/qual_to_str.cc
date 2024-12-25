#include "art_modern_dtypes.hh"
#include "utils/seq_utils.hh"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace labw::art_modern;

void bench(const int run_times, const int rlen)
{
    std::cout << "run_times: " << run_times << " rlen: " << rlen << std::endl;
    std::vector<am_qual_t> q;
    q.reserve(rlen);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 40);
    for (int i = 0; i < rlen; i++) {
        q.emplace_back(dis(gen));
    }
    std::chrono::high_resolution_clock::time_point start, end;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_avx2(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "AVX2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_sse2(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "SSE2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_mmx(q.data(), rlen);
    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "MMX: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Scala (for loop): " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_foreach(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Scala (std::for_each): " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << std::endl;

    if (qual_to_str_sse2(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        throw std::exception();
    }

    if (qual_to_str_mmx(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        throw std::exception();
    }

    if (qual_to_str_for_loop(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        throw std::exception();
    }

    if (qual_to_str_avx2(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        throw std::exception();
    }
}

int main()
{
    bench(1000000, 36);
    bench(1000000, 100);
    bench(1000000, 150);
    bench(1000000, 200);
    bench(1000000, 300);
    bench(100000, 1000);
    bench(10000, 10000);
    return 0;
}
