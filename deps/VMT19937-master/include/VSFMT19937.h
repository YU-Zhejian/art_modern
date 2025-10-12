#pragma once

#include "SIMD.h"
#include "jump_matrix.h"
#include "Params.h"

#include <cstdint>
#include <cstddef>

namespace Details {

template <size_t RegisterBitLen, size_t RegisterBitLenImpl>
class VSFMT19937Base
{
    static const size_t s_nBits = SFMT19937Params::s_nBits;
    static const size_t s_wordSizeBits = SFMT19937Params::s_wordSizeBits;

    static_assert(RegisterBitLen >= s_wordSizeBits);

public:
    static const int s_N = SFMT19937Params::s_N;     // 156
    static_assert(s_N == 156);

    static const size_t s_regLenBits = RegisterBitLen;
    static const size_t s_regLenImplBits = RegisterBitLenImpl;
    static const size_t s_nStates = RegisterBitLen / s_wordSizeBits;
    static const size_t s_n32inReg = RegisterBitLen / 32;
    static const size_t s_n32InOneWord = s_wordSizeBits / 32;            // 4
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;          // 624
    const static size_t s_n32InFullState = s_n32InOneState * s_nStates;  // 624 * nStates
    const static size_t s_nMatrixBits = SFMT19937Params::s_nMatrixBits;            // 19968

    typedef BinaryMatrix<s_nMatrixBits> matrix_t;

private:
    const static size_t s_regLenWords = s_regLenBits / s_wordSizeBits;  // FIXME: review this definition

    typedef SimdRegister<s_regLenBits, RegisterBitLenImpl> XV;

protected:
    alignas(64) uint32_t m_state[s_n32InFullState];    // the array of state vectors

private:
    // This data members is necessary only if QueryMode!=QM_StateSize
    const uint32_t* const m_state_end;
    const uint32_t* m_prnd;

    template <typename XVCst>
    static FORCE_INLINE XV advance1(const XV& xA, const XV& xB, const XV& xC, const XV& xD, const XVCst& bMask)
    {
        const uint32_t SFMT_SL1 = 18;
        const uint32_t SFMT_SL2 = 1;
        const uint32_t SFMT_SR1 = 11;
        const uint32_t SFMT_SR2 = 1;

        XV y(xB >> SFMT_SR1);
        XV z(XV::template shr128<SFMT_SR2>(xC));
        XV v(xD << SFMT_SL1);
        z = z ^ xA;
        z = z ^ v;
        XV x(XV::template shl128<SFMT_SL2>(xA));
        y = y & bMask;
        z = z ^ x;
        return (z ^ y);
    }

    template <int nIter, int JA, int JB, typename XVCst>
    static FORCE_INLINE void unroll(uint32_t* p, XV& xC, XV& xD, const XVCst& bMask)
    {
        if constexpr (nIter > 0) {
            XV xA(p + JA * s_n32inReg);
            XV xB(p + JB * s_n32inReg);
            XV tmp = advance1(xA, xB, xC, xD, bMask);
            xC = xD;
            xD = tmp;
            tmp.template store<true>(p + JA * s_n32inReg);
            unroll<nIter - 1, JA + 1, JB + 1>(p, xC, xD, bMask);
        }
    }

    template <int nUnroll, int nIter, int JB, typename XVCst>
    static FORCE_INLINE void advanceLoop(uint32_t*& p, XV& xC, XV& xD, const XVCst& bMask)
    {
        // unroll main iterations
        const size_t nMainIter = nIter / nUnroll;
        if constexpr (nMainIter) {
            auto pend = p + nMainIter * nUnroll * s_n32inReg;
            // unroll the loop in blocks of UnrollBlkSize
            do {
                unroll<nUnroll, 0, JB>(p, xC, xD, bMask);
                p += nUnroll * s_n32inReg;
            } while (p != pend);
        }

        // unroll residual iterations (if any)
        const size_t nResIter = nIter % nUnroll;
        if constexpr (nResIter) {
            unroll<nResIter, 0, JB>(p, xC, xD, bMask);
            p += nResIter * s_n32inReg;
        }
    }

    NO_INLINE void refill()
    {
        // Create local copy of the constants and pass them to the function as arguments.
        // Since all functions invoked from here are forced inline, the function arguments
        // will not be passed as arguments via the stack, but reside in CPU registers
        SimdRegister<std::max<size_t>(128, s_regLenImplBits), s_regLenImplBits>
            bMask(SFMT19937Params::s_SFMT_MSK1, SFMT19937Params::s_SFMT_MSK2, SFMT19937Params::s_SFMT_MSK3, SFMT19937Params::s_SFMT_MSK4);

        const int s_M = SFMT19937Params::s_M;

        // local variables
        uint32_t* stCur = m_state;
        XV xC(stCur + (s_N - 2) * s_n32inReg);
        XV xD(stCur + (s_N - 1) * s_n32inReg);

        // unroll first part of the loop: (N-M) iterations
        advanceLoop<2, s_N - s_M, s_M>(stCur, xC, xD, bMask);

        // unroll second part of the loop: M iterations
        advanceLoop<2, s_M, s_M - s_N>(stCur, xC, xD, bMask);

        m_prnd = begin();
    }

    const uint32_t* begin() const
    {
        return m_state;
    }

    const uint32_t* end() const
    {
        return m_state + s_n32InFullState;
    }

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    // return the absolute index in the state vector
    // of the 32-bit word `w32RelIndex` belonging to state `stateIndex`
    static size_t w32AbsIndex(size_t w32RelIndex, size_t stateIndex = 0)
    {
        size_t w128RelIndex = w32RelIndex / s_n32InOneWord;
        size_t w128RelOffset = w32RelIndex % s_n32InOneWord;
        return w128RelIndex * s_n32inReg + s_n32InOneWord * stateIndex + w128RelOffset;
    }

    // This function ensures the period of 2^{MEXP}
    void ensure_period()
    {
        const uint32_t SFMT_PARITY1 = 0x00000001U;
        const uint32_t SFMT_PARITY2 = 0x00000000U;
        const uint32_t SFMT_PARITY3 = 0x00000000U;
        const uint32_t SFMT_PARITY4 = 0x13c9e684U;

        const uint32_t parity[4] = { SFMT_PARITY1, SFMT_PARITY2, SFMT_PARITY3, SFMT_PARITY4 };

        uint32_t inner = 0;
        uint32_t* psfmt32 = m_state;

        for (size_t i = 0; i < 4; i++)
            inner ^= psfmt32[i] & parity[i];

        for (size_t i = 16; i > 0; i >>= 1)
            inner ^= inner >> i;

        inner &= 1;

        // check if OK
        if (inner == 1)
            return;

        // check NG, and modification
        for (size_t i = 0; i < 4; i++) {
            uint32_t work = 1;
            for (size_t j = 0; j < 32; j++) {
                if ((work & parity[i]) != 0) {
                    psfmt32[i] ^= work;
                    return;
                }
                work = work << 1;
            }
        }
    }

    // convenience used in the initialization
    static uint32_t func1(uint32_t x)
    {
        return (x ^ (x >> 27)) * (uint32_t)1664525UL;
    }

    // convenience used in the initialization
    static uint32_t func2(uint32_t x)
    {
        return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
    }

    void getBlock16(uint32_t* dst)
    {
        std::copy(m_prnd, m_prnd + 16, dst);
        m_prnd += 16;
    }

protected:

    // extract one of the interleaved state vectors and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        cubeToMatrix<s_nStates, s_n32InOneWord, uint32_t>(pdst, m_state, s_N, stateIndex);
    }

    // store vector into the interleaved elements of the state vector
    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        matrixToCube<s_nStates, s_n32InOneWord, uint32_t>(m_state, psrc, s_N, stateIndex);
    }

    // Initializes the internal state array with a 32-bit integer seed.
    void reinitMainState(uint32_t seed)
    {
        uint32_t* psfmt32 = begin();

        psfmt32[0] = seed;
        for (size_t i = 1; i < s_n32InOneState; i++) {
            psfmt32[w32AbsIndex(i)] = 1812433253UL * (psfmt32[w32AbsIndex(i - 1)]
                ^ (psfmt32[w32AbsIndex(i - 1)] >> 30))
                + i;
        }
        ensure_period();
    }

    // Initializes the internal state array, with an array of 32-bit integers used as the seeds
    void reinitMainState(const uint32_t* init_key, uint32_t key_length)
    {
        uint32_t i, j;
        uint32_t* psfmt32 = m_state;

        const size_t lag = 11;
        const size_t mid = (s_n32InOneState - lag) / 2;

        for (size_t w = 0; w < s_N; ++w)
            for (size_t j = 0; j < s_n32InOneWord; ++j)
                psfmt32[w * s_n32inReg + j] = uint32_t(0x8b8b8b8b);
        size_t count = (key_length + 1 > s_n32InOneState) ? key_length + 1 : s_n32InOneState;
        uint32_t r = func1(psfmt32[w32AbsIndex(0)] ^ psfmt32[w32AbsIndex(mid)] ^ psfmt32[w32AbsIndex(s_n32InOneState - 1)]);
        psfmt32[w32AbsIndex(mid)] += r;
        r += key_length;
        psfmt32[w32AbsIndex(mid + lag)] += r;
        psfmt32[w32AbsIndex(0)] = r;

        count--;
        for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
            r = func1(psfmt32[w32AbsIndex(i)] ^ psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^ psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] += r;
            r += init_key[j] + i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] += r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }
        for (; j < count; j++) {
            r = func1(psfmt32[w32AbsIndex(i)] ^ psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^ psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] += r;
            r += i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] += r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }
        for (j = 0; j < s_n32InOneState; j++) {
            r = func2(psfmt32[w32AbsIndex(i)] + psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] + psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^= r;
            r -= i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] ^= r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }

        ensure_period();
    }

    void reinitPointers()
    {
        m_prnd = end();
    }

    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        if (m_prnd != end())
            return *m_prnd++;

        refill();

        return *m_prnd++;
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        if (m_prnd != end()) {
            getBlock16(dst);
        }
        else {
            refill();
            getBlock16(dst);
        }
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        refill();
        std::copy(begin(), end(), dst);
    }

    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        if (size_t nAvail = std::distance(m_prnd, end());  nAvail <= n) {
            std::copy_n(m_prnd, nAvail, dst);
            n -= nAvail;
            dst += nAvail;
            refill();
        }
        else {
            std::copy_n(m_prnd, n, dst);
            m_prnd += n;
            return;
        }

        while (n >= s_n32InFullState) {
            std::copy_n(m_prnd, s_n32InFullState, dst);
            n -= s_n32InFullState;
            dst += s_n32InFullState;
            refill();
        }

        std::copy_n(begin(), n, dst);
        m_prnd += n;
    }

public:

    // constructors
    VSFMT19937Base()
        : m_state_end(m_state + s_n32InFullState)
        , m_prnd(nullptr)
    {}

}; // VSFMT19937Base

} // namespace Details

