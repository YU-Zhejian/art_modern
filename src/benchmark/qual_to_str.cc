#include "art_modern_constants.hh"
#include "art_modern_dtypes.hh"
#include <chrono>
#include <iostream>
#include <random>
#include <simde/simde-common.h>
#include <simde/x86/sse2.h>
#include <string>
#include <vector>

using namespace labw::art_modern;

std::string qual_to_str_sse2(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    size_t i = 0;
#ifdef SIMDE_X86_SSE2_NATIVE
    const size_t num_elements_per_simd = 16; // SSE2 processes 16 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    simde__m128i phred_offset_vec = simde_mm_set1_epi8(static_cast<uint8_t>(PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        simde__m128i qual_vec = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(&qual[i]));
        simde__m128i result_vec = simde_mm_add_epi8(qual_vec, phred_offset_vec);
        simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(&retq[i]), result_vec);
    }
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

std::string qual_to_str(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    // Handle the remaining elements that do not fit into a full SIMD register
    for (size_t i = 0; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

int main()
{
    const int run_times = 100000;
    const int rlen = 151;
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
        volatile auto s = qual_to_str_sse2(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "SSE2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str(q.data(), rlen);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Scala: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    if (qual_to_str_sse2(q.data(), rlen) != qual_to_str(q.data(), rlen)) {
        throw std::exception();
    }
    return 0;
}