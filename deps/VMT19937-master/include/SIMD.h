#pragma once

#include "simd_config.h"
#include "portable.h"

#include <algorithm>

namespace Details {

template <size_t NumBits, size_t NumBitsImpl>
struct SimdRegister
{
    static const size_t s_nRegBits = NumBits;
    static const size_t s_nImplBits = NumBitsImpl;
    static_assert(NumBits >= 2 * NumBitsImpl);
private:
    static_assert((NumBits / NumBitsImpl > 1) && (NumBits % NumBitsImpl == 0), "NumBits must be a multiple of NumBitsImpl");
    static const size_t M = NumBits / NumBitsImpl;
    static const size_t N32 = NumBits / 32;
    static const size_t N128 = NumBits / 128;
    typedef SimdRegister<NumBitsImpl, NumBitsImpl> XVImpl;

    struct Aux
    {
        Aux() : ar{ {} } {}
        Aux(XVImpl v) { std::fill_n(ar, M, v); }
        Aux(uint32_t v) : Aux(XVImpl(v)) {}
        XVImpl& operator[](size_t i) { return ar[i]; }
        const XVImpl& operator[](size_t i) const { return ar[i]; }
        const XVImpl* begin() const { return ar; }
        XVImpl* begin() { return ar; }
        const XVImpl* end() const { return ar + M; }
        XVImpl* end() { return ar + M; }
    private:
        XVImpl ar[M];
    };

public:
    Aux m_v;

    typedef SimdRegister<NumBits, NumBitsImpl> XV;

    SimdRegister() {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(v) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        if constexpr (NumBitsImpl == 32) {
            for (size_t i = 0; i < N128; ++i) {
                m_v[0 + 4 * i] = v0;
                m_v[1 + 4 * i] = v1;
                m_v[2 + 4 * i] = v2;
                m_v[3 + 4 * i] = v3;
            }
        }
        else {
            m_v = Aux(XVImpl(v0, v1, v2, v3));
        }
    }
    FORCE_INLINE SimdRegister(const uint32_t* p)
    {
        for (size_t i = 0; i < M; ++i, p += sizeof(XVImpl) / sizeof(uint32_t))
            m_v[i] = XVImpl(p);
    }
    FORCE_INLINE SimdRegister(const Aux& v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* p)
    {
        for (size_t i = 0; i < M; ++i, p += sizeof(XVImpl) / sizeof(uint32_t))
            m_v[i].template store<A>(p);
    }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v[i] = a.m_v[i] & b.m_v[i]; return r; }
    template <typename XVI>
    friend FORCE_INLINE XV operator&(const XV& a, const XVI& b)
    {
        XV r;
        if constexpr (XVI::s_nRegBits > XVI::s_nImplBits) {
            const size_t N = XVI::s_nRegBits / XVI::s_nImplBits;
            for (size_t i = 0; i < M / N; ++i)
                for (size_t j = 0; j < N; ++j)
                    r.m_v[i * N + j] = a.m_v[i * N + j] & b.m_v[j];
        }
        else {
            for (size_t i = 0; i < M; ++i)
                r.m_v[i] = a.m_v[i] & b;
        }
        return r;
    }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v[i] = a.m_v[i] ^ b.m_v[i]; return r; }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v[i] = a.m_v[i] | b.m_v[i]; return r; }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { XV r; for (size_t i = 0; i < M; ++i) r.m_v[i] = a.m_v[i] << n; return r; }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { XV r; for (size_t i = 0; i < M; ++i) r.m_v[i] = a.m_v[i] >> n; return r; }

    template <int nBytes>
    FORCE_INLINE static XV shl128(const XV& a)
    {
        static_assert(N128 > 0);
        XV r;
        if constexpr (NumBitsImpl == 32) {
            for (size_t s = 0; s < N128; ++s) {
                XVImpl current = a.m_v[4 * s];
                r.m_v[4 * s] = current << 8;
                for (size_t j = 1; j < 4; ++j) {
                    XVImpl prev = current;
                    current = a.m_v[4 * s + j];
                    r.m_v[4 * s + j] = (current << 8) | (prev >> 24);
                }
            }
        }
        else {
            for (size_t i = 0; i < M; ++i)
                r.m_v[i] = XVImpl::template shl128<nBytes>(a.m_v[i]);
        }
        return r;
    }

    template <int nBytes>
    FORCE_INLINE static XV shr128(const XV& a)
    {
        static_assert(N128 > 0);
        XV r;
        if constexpr (NumBitsImpl == 32) {
            for (size_t s = 0; s < N128; ++s) {
                XVImpl current = a.m_v[4 * s];
                for (size_t j = 0; j < 3; ++j) {
                    XVImpl next = a.m_v[4 * s +  j + 1];
                    r.m_v[4 * s + j] = (current >> 8) | (next << 24);
                    current = next;
                }
                r.m_v[4 * s + 3] = current >> 8;
            }
        }
        else {
            for (size_t i = 0; i < M; ++i)
                r.m_v[i] = XVImpl::template shr128<nBytes>(a.m_v[i]);
        }
        return r;
    }

    FORCE_INLINE void broadcastLo128()
    {
        static_assert(N128 > 0);
        if constexpr (NumBitsImpl == 32) {
            for (size_t i = 0; i < N128; ++i)
                for (size_t j = 0; j < 4; ++j)
                    m_v[4 * i + j] = m_v[j];
        }
        else {
            for (size_t i = 1; i < M; ++i)
                m_v[i] = m_v[0];
        }
    }

    bool eq(const XV& rhs) const
    {
        for (size_t i = 0; i < M; ++i)
            if (m_v[i] != rhs.m_v[i])
                return false;
        return true;
    }

    FORCE_INLINE static XV zero()
    {
        XVImpl z(0);
        XV r;
        for (auto& v : r.m_v)
            v = z;
        return r;
    }

    template <typename XVI>
    FORCE_INLINE XV ifOddCst32ElseZero(const XVI& value) const
    {
        XV r;
        for (size_t i = 0; i < M; ++i)
            r.m_v[i] = m_v[i].ifOddCst32ElseZero(value);
        return r;
    }

    FORCE_INLINE uint8_t parity() const
    {
        XVImpl temp(m_v[0]);
        for (size_t i = 1; i < M; ++i)
            temp = temp ^ m_v[i];
        return temp.parity();
    }
};


template <>
struct SimdRegister<32, 32>
{
    static const size_t s_nRegBits = 32;
    static const size_t s_nImplBits = 32;

    uint32_t m_v;

    typedef SimdRegister<32, 32> XV;

    SimdRegister() : m_v(0) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(*(const uint32_t*)p) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { *dst = m_v; }

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return a.m_v > b.m_v ? uint32_t(-1) : uint32_t(0); }
    friend FORCE_INLINE XV operator>>(const XV a, int n) { return uint32_t(a.m_v >> n); }
    friend FORCE_INLINE XV operator<<(const XV a, int n) { return uint32_t(a.m_v << n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
#if 0 && defined(__GNUC__) && defined(__x86_64__)
        // force the use of cmov with gcc
        uint32_t z;
        __asm__(
            "mov %[a], %[z]\n"
            "and $0x1, %[z]\n"
            "cmovne %[b], %[z]\n"
            : [z] "=r"(z)
            : [a] "r"(m_v), [b] "r"(cst32.m_v)
            : "cc"
        );
        return z;
#elif defined(_MSC_VER)
        const uint32_t lowestBit = m_v & 0x1;
        return lowestBit ? cst32 : zero();
#else
        const uint32_t x[2] = { 0, cst32.m_v };
        const uint32_t lowestBit = m_v & 0x1;
        return x[lowestBit];
#endif
    }

    static FORCE_INLINE XV zero() { return uint32_t(0); }

    uint8_t parity() const { return popcnt(m_v) % 2; }
};


/*
template <>
struct MAY_ALIAS SimdRegister<64>
{
    uint64_t m_v;

    typedef SimdRegister<64> XV;

    SimdRegister() : m_v(0) {}
    SimdRegister(const void* p) : m_v(*(const uint64_t*)p) {}
    SimdRegister(uint32_t v) : m_v(v | (uint64_t(v) << 32)) {}
    SimdRegister(uint64_t v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return a.m_v | b.m_v; }
    friend FORCE_INLINE XV operator>>(const XV& a, int n)
    {
        uint64_t mask = (uint64_t(0xFFFFFFFF) << 32) | ((uint32_t(1) << (32 - n)) - 1);
        return uint64_t(a.m_v >> n) & mask;
    }
    friend FORCE_INLINE XV operator<<(const XV& a, int n)
    {
        //uint64_t mask = uint64_t(-1) & ~(((uint64_t(1) << n) - 1) << 32);
        return uint64_t(a.m_v << n);
    }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV& value) const
    {
        const uint64_t maskHi = m_v & (uint64_t(1) << 32) ? (uint64_t(0xFFFFFFFF) << 32) : 0;
        const uint64_t maskLo = m_v & 0x1 ? uint64_t(0xFFFFFFFF) : 0;
        return value.m_v & (maskHi | maskLo);
    }

    static FORCE_INLINE XV zero() { return uint64_t(0); }

    uint8_t parity() const { return popcnt(m_v) % 2; }
};
*/

#if SIMD_N_BITS>=128
template <>
struct SimdRegister<128, 128>
{
    static const size_t s_nRegBits = 128;
    static const size_t s_nImplBits = 128;

    __m128i m_v;

    typedef SimdRegister<128, 128> XV;

    SimdRegister() : m_v(_mm_undefined_si128()) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm_setr_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm_load_si128((const __m128i*)p)) {}
    FORCE_INLINE SimdRegister(__m128i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm_store_si128((__m128i*)dst, m_v); else _mm_storeu_si128((__m128i*)dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return _mm_test_all_ones(_mm_cmpeq_epi32(m_v, rhs.m_v)); }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm_bslli_si128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm_bsrli_si128(a.m_v, n); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV& cst32) const
    {
#if 0
        const __m128 z = _mm_setzero_ps();
        const __m128 lowestBit = _mm_castsi128_ps(_mm_slli_epi32(m_v, 31));
        return _mm_castps_si128(_mm_blendv_ps(z, _mm_castsi128_ps(cst32.m_v), lowestBit));
#else
        const __m128i z = zero().m_v;
        const __m128i lowestBit = _mm_slli_epi32(m_v, 31);
        const __m128i isOdd = _mm_cmpgt_epi32(z, lowestBit);
        return _mm_and_si128(isOdd, cst32.m_v);
#endif
    }

    uint8_t parity() const
    {
        __m128i hi(_mm_shuffle_epi32(m_v, 2 | (3 << 2)));
        __m128i mix = _mm_xor_si128(hi, m_v);
        uint64_t d = _mm_extract_epi64(mix, 0);
        return popcnt(d) & 1;
    }
};
#endif

#if SIMD_N_BITS>=256
template <>
struct SimdRegister<256, 256>
{
    static const size_t s_nRegBits = 256;
    static const size_t s_nImplBits = 256;

    __m256i m_v;

    typedef SimdRegister<256, 256> XV;

    SimdRegister() : m_v(_mm256_undefined_si256()) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm256_load_si256((const __m256i*)p)) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm256_set1_epi32(v)) {}
    FORCE_INLINE FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        __m128i tmp = _mm_setr_epi32(v0, v1, v2, v3);
        m_v = _mm256_set_m128i(tmp, tmp);
    }
    SimdRegister(__m256i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm256_store_si256((__m256i*)dst, m_v); else _mm256_storeu_si256((__m256i*)dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm256_or_si256(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm256_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm256_srli_epi32(a.m_v, n); }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm256_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm256_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm256_broadcastsi128_si256(_mm256_castsi256_si128(m_v)); }

    static FORCE_INLINE XV zero() { return _mm256_setzero_si256(); }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV& cst32) const
    {
        const __m256i z = zero().m_v;
        const __m256i lowestBit = _mm256_slli_epi32(m_v, 31); // move least significant bit to most significant bit
#if 0
        const __m256i isOdd = _mm256_cmpgt_epi32(z, lowestBit);
        return _mm256_and_si256(value.m_v, isOdd);
#else
        const __m256 mask = _mm256_castsi256_ps(lowestBit);
        return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(z), _mm256_castsi256_ps(cst32.m_v), mask));
#endif
    }

    uint8_t parity() const
    {
        __m128i hi(_mm256_extracti128_si256(m_v, 1));
        __m128i lo(_mm256_castsi256_si128(m_v));
        return SimdRegister<128, 128>(_mm_xor_si128(lo, hi)).parity();
    }
};
#endif

#if SIMD_N_BITS>=512
template <>
struct SimdRegister<512, 512>
{
    static const size_t s_nRegBits = 512;
    static const size_t s_nImplBits = 512;

    __m512i m_v;

    typedef SimdRegister<512, 512> XV;

    SimdRegister() : m_v(_mm512_undefined_si512()) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm512_load_si512((const __m512i*)p)) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm512_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm512_setr4_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(__m512i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm512_store_si512((__m512i*)dst, m_v); else _mm512_storeu_si512((__m512i*)dst, m_v); }

    //    static FORCE_INLINE XV load(const void* p) { return _mm512_loadu_si512((const __m512i*) p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm512_or_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm512_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm512_srli_epi32(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm512_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm512_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm512_broadcast_i32x4(_mm512_castsi512_si128(m_v)); }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV& cst32) const
    {
        const __m512i lowestBit = _mm512_slli_epi32(m_v, 31); // move least significant bit to most significant bit
        const __mmask16 isOdd = _mm512_movepi32_mask(lowestBit);
        return _mm512_maskz_mov_epi32(isOdd, cst32.m_v);
    }

    static FORCE_INLINE XV zero() { return _mm512_setzero_si512(); }

    uint8_t parity() const
    {
        __m256i hi(_mm512_extracti64x4_epi64(m_v, 1));
        __m256i lo(_mm512_castsi512_si256(m_v));
        return SimdRegister<256, 256>(_mm256_xor_si256(lo, hi)).parity();
    }
};
#endif

} // namespace Details

