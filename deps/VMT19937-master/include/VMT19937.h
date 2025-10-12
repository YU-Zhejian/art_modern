#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>

namespace Details {

template <size_t RegisterBitLen, size_t RegisterBitLenImpl>
class VMT19937Base
{

    static const size_t s_nBits = MT19937Params::s_nBits;
    static const size_t s_wordSizeBits = MT19937Params::s_wordSizeBits;

    static_assert(RegisterBitLen >= s_wordSizeBits);

public:

    static const int s_N = MT19937Params::s_N;     // 624
    static_assert(s_N == 624);

    static const size_t s_regLenBits = RegisterBitLen;
    static const size_t s_regLenImplBits = RegisterBitLenImpl;
    static const size_t s_nStates = RegisterBitLen / s_wordSizeBits;
    static const size_t s_n32inReg = RegisterBitLen / 32;
    static const size_t s_n32InOneWord = s_wordSizeBits / 32;            // 1
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;          // 624
    const static size_t s_n32InFullState = s_n32InOneState * s_nStates;  // 624 * nStates
    const static size_t s_nMatrixBits = MT19937Params::s_nMatrixBits;

    typedef MT19937Matrix matrix_t;

private:
    const static size_t s_regLenWords = s_regLenBits / s_wordSizeBits;  // FIXME: review this definition
    typedef SimdRegister<s_regLenBits, RegisterBitLenImpl> XV;

    static const uint32_t s_cacheLineBits = 64*8;
    static const uint32_t s_n32InRndCache = std::max<uint32_t>(s_cacheLineBits, 4*RegisterBitLenImpl) / 32;
    static_assert(s_n32InFullState % s_n32InRndCache == 0, "full state size not divisible by cache size");
    static_assert(s_n32InRndCache % (64/4) == 0, "cache size is not divisible by cache line size");

protected:
    alignas(64) uint32_t m_state[s_N * s_n32inReg];    // the array of state vectors

private:
    // This data members is necessary only if QueryMode==QM_Scalar
    alignas(64) uint32_t m_rnd[s_n32InRndCache]; // buffer of tempered numbers
    const uint32_t* m_prnd, * const m_prndEnd;

    // This data members are redundant if QueryMode==QM_StateSize
    const uint32_t*m_pst, * const m_pstEnd;    // m_pos==m_pstEnd means the state vector has been consumed and need to be regenerated

private:

    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask1(MT19937Params::s_temperMask1), m_mask2(MT19937Params::s_temperMask2) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    struct RefillCst
    {
        typedef SimdRegister<RegisterBitLenImpl, RegisterBitLenImpl> XVI;
        RefillCst() : m_upperMask(MT19937Params::s_upperMask), m_lowerMask(MT19937Params::s_lowerMask), m_matrixA(MT19937Params::s_matrixA) {}
        const XVI m_upperMask;
        const XVI m_lowerMask;
        const XVI m_matrixA;
    };

    template <typename XVI, typename M>
    static FORCE_INLINE XVI temper(XVI y, const M& masks)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & masks.m_mask1);
        y = y ^ ((y << 15) & masks.m_mask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <bool Aligned>
    static FORCE_INLINE void temperRefillBlock_(const uint32_t * __restrict st, uint32_t * __restrict dst)
    {
        const size_t LR = std::min<size_t>(s_regLenImplBits * 4, s_n32InRndCache * 32);
        typedef SimdRegister<LR, s_regLenImplBits> XVmax;
        const size_t n32PerIteration = LR / 32;
        static_assert(n32PerIteration <= s_n32InRndCache);
        static_assert(s_n32InRndCache % n32PerIteration == 0);
        const size_t nIterations = s_n32InRndCache / n32PerIteration;

        TemperCst<XVmax> cst{};

        for (size_t i = 0; i < nIterations; ++i) {
            XVmax tmp = temper(XVmax(st), cst);
            tmp.template store<Aligned>(dst);
            dst += n32PerIteration;
            st += n32PerIteration;
        }
    }

    template <bool Aligned>
    static FORCE_INLINE void temperRefillBlock(const uint32_t*& st_, uint32_t* dst)
    {
        temperRefillBlock_<Aligned>(st_, dst);
        st_ += s_n32InRndCache;
    }

    static FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm, const RefillCst& masks)
    {
        XV y = (s & masks.m_upperMask) | (sp & masks.m_lowerMask);
        // y and sp are either both even or both odd,
        // hence in the next line we can check if sp is odd
        // so that the operation is independent on the calculation of y
        // and the compiler is free to rearrange the code
        XV r = sm ^ (y >> 1) ^ sp.ifOddCst32ElseZero(masks.m_matrixA);
        return r;
    }

    //template <int nIter>
    static FORCE_INLINE void iteration(uint32_t* p, XV& x0, int J0, int J1, int JM, const RefillCst& masks)
    {
        XV x1(p + J1 * s_n32inReg);
        XV xM(p + JM * s_n32inReg);
        XV tmp = advance1(x0, x1, xM, masks);
        tmp.template store<true>(p + J0 * s_n32inReg);
        x0 = x1;
    }

    template <int...Is>
    static FORCE_INLINE uint32_t* advanceLoop(size_t nBlkIter, uint32_t* p, XV& x0, int J1, int JM, const RefillCst& masks, std::integer_sequence<int, Is...>&&)
    {
        const size_t nIterPerBlk = sizeof...(Is);
        if constexpr (nIterPerBlk) {
            const size_t n32PerBlk = nIterPerBlk * s_n32inReg;
            auto pend = p + nBlkIter * n32PerBlk;
            do {
                (iteration(p, x0, 0 + Is, J1 + Is, JM + Is, masks), ...);
                p += n32PerBlk;
            } while (p != pend);
            return pend;
        }
        else
            return p;
    }

    void NO_INLINE refill()
    {
        uint32_t* stCur = m_state;

        const int N = s_N;
        const int M = Details::MT19937Params::s_M;
        static_assert(N == 624 && M == 397, "unrolling designed for these parameters");

        // Create local copy of the constants and pass them to the function as arguments.
        // Since all functions invoked from here are forced inline, the function arguments
        // will not be passed as arguments via the stack, but reside in CPU registers
        const RefillCst masks;  // use default constructor

        XV x0(stCur);

        const size_t nUnroll = 2;

        // unroll first part of the loop (N-M) iterations
        const size_t n1 = (N - M) / nUnroll;
        const size_t r1 = (N - M) % nUnroll;
        stCur = advanceLoop(n1, stCur, x0, 1, M, masks, std::make_integer_sequence<int, nUnroll>{});
        stCur = advanceLoop(r1, stCur, x0, 1, M, masks, std::make_integer_sequence<int, r1>{});

        // unroll second part of the loop (M-1) iterations
        const size_t n2 = (M - 1) / nUnroll;
        const size_t r2 = (M - 1) % nUnroll;
        stCur = advanceLoop(n2, stCur, x0, 1, M - N, masks, std::make_integer_sequence<int, nUnroll>{});
        stCur = advanceLoop(r2, stCur, x0, 1, M - N, masks, std::make_integer_sequence<int, r2>{});

        // last iteration
        advanceLoop(1, stCur, x0, 1 - N, M - N, masks, std::make_integer_sequence<int, 1>{});

        m_pst = begin();
    }

    const uint32_t* begin() const
    {
        return m_state;
    }

    const uint32_t* end() const
    {
        return m_pstEnd;
    }

    const uint32_t* beginRnd() const
    {
        return m_rnd;
    }

    const uint32_t* endRnd() const
    {
        return m_prndEnd;
    }


    uint32_t& scalarState(uint32_t scalarIndex)
    {
        return m_state[scalarIndex * s_regLenWords];
    }

    uint32_t scalarState(uint32_t scalarIndex) const
    {
        return m_state[scalarIndex * s_regLenWords];
    }

    // initializes the first state with a seed
    void __reinit(uint32_t s)
    {
        const uint32_t mask = uint32_t(1812433253UL);
        uint32_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < s_N; i++)
            prev = scalarState(i) = (mask * (prev ^ (prev >> 30)) + i);
    }

protected:
    void reinitPointers()
    {
        m_pst = m_pstEnd;
        m_prnd = (const uint32_t*)(((uint8_t*)m_rnd) + sizeof(m_rnd));
    }

    // extract one of the interleaved state vectors, shift it left by 31 bits and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint32_t* pstate = m_state;
        pdst[0] = pstate[stateIndex] >> 31;
        for (size_t i = 1; i < s_N; ++i) {
            uint32_t word = pstate[i * s_regLenWords + stateIndex];
            pdst[i - 1] |= word << 1;
            pdst[i] = word >> 31;
        }
    }

    // shift vector psr to the right by 31 bit and store into the interleaved elements of the state vector
    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        uint32_t* pstate = m_state;
        const uint32_t* pw = (const uint32_t*)psrc;
        pstate[stateIndex] = 0;
        size_t w;
        for (w = 0; w < s_N - 1; ++w) {
            uint32_t word = pw[w];
            pstate[w * s_regLenWords + stateIndex] |= word << 31;
            pstate[(w + 1) * s_regLenWords + stateIndex] = word >> 1;
        }
        pstate[w * s_regLenWords + stateIndex] |= pw[w] << 31;
    }

    // initializes m_state[s_N] with a seed
    void reinitMainState(uint32_t s)
    {
        __reinit(s);
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    void reinitMainState(const uint32_t* seeds, uint32_t nSeeds)
    {
        __reinit(uint32_t(19650218));
        uint32_t i = 1, j = 0;
        uint32_t k = (s_N > nSeeds ? s_N : nSeeds);
        for (; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * uint32_t(1664525)))
                + seeds[j] + j; // non linear
            //m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++; j++;
            if (i >= s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
        }
        for (k = s_N - 1; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * uint32_t(1566083941)))
                - i; // non linear 
            //m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++;
            if (i >= s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
        }

        scalarState(0) = uint32_t(0x80000000); // MSB is 1; assuring non-zero initial array
    }

    // generates a random number on [0,0xffffffff] interval
    FORCE_INLINE uint32_t genrand_uint32()
    {
        if (m_prnd != endRnd())
            return *m_prnd++;

        if (m_pst != m_pstEnd)
            /* do nothing*/; // most likely case first
        else {
            refill();
        }

        temperRefillBlock<true>(m_pst, m_rnd);
        m_prnd = beginRnd() + 1;

        return m_rnd[0];
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_uint32_blk16(uint32_t* dst)
    {
        if constexpr (s_n32InRndCache > 16) {
            if (m_prnd != endRnd()) {
                const uint32_t* e = m_prnd + 16;
                std::copy(m_prnd, e, dst);
                m_prnd = e;
                return;
            }
        }

        if (m_pst != m_pstEnd)
            /* do nothing*/; // most likely case first
        else
            refill();

        if constexpr (s_n32InRndCache > 16) {
            temperRefillBlock<true>(m_pst, m_rnd);
            uint32_t* e = m_rnd + 16;
            std::copy(m_rnd, e, dst);
            m_prnd = e;
        }
        else {
            temperRefillBlock<false>(m_pst, dst);
        }
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        refill();
        const uint32_t* pst = m_state;
        for (size_t i = 0; i < s_n32InFullState / s_n32InRndCache; ++i, dst += s_n32InRndCache)
            temperRefillBlock<false>(pst, dst);
    }

    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        if (size_t nAvailInRnd = std::distance<const uint32_t*>(m_prnd, endRnd()); nAvailInRnd < n) {
            std::copy_n(m_prnd, nAvailInRnd, dst);
            n -= nAvailInRnd;
            dst += nAvailInRnd;
        }
        else {
            std::copy_n(m_prnd, n, dst);
            m_prnd += n;
            return;
        }

        // we are now aligned with rnd block
        if (size_t nAvailInState = std::distance(m_pst, end()); nAvailInState < n) {
            if (nAvailInState) {
                n -= nAvailInState;
                do {
                    temperRefillBlock<false>(m_pst, dst);
                    dst += s_n32InRndCache;
                    nAvailInState -= s_n32InRndCache;
                } while (nAvailInState);
            }
        }
        else {
            size_t nFullRndBlocks = n / s_n32InRndCache;
            for (size_t i = 0; i < nFullRndBlocks; ++i) {
                temperRefillBlock<false>(m_pst, dst);
                dst += s_n32InRndCache;
            }
            n = n % s_n32InRndCache;
            if (n) {
                temperRefillBlock<true>(m_pst, m_rnd);
                std::copy_n(m_rnd, n, dst);
                m_prnd = m_rnd + n;
            }
            else {
                m_prnd = endRnd();
            }
            return;
        }

        // we are now aligned with the state vector

        size_t nFullStates = n / s_n32InFullState;
        while (nFullStates--) {
            genrand_uint32_stateBlk(dst);
            dst += s_n32InFullState;
            n -= s_n32InFullState;
        }
        m_pst = end();

        if (n > 0) {
            refill();

            while (n > s_n32InRndCache) {
                temperRefillBlock<false>(m_pst, dst);
                n -= s_n32InRndCache;
                dst += s_n32InRndCache;
            }

            if (n) {
                temperRefillBlock<true>(m_pst, m_rnd);
                std::copy_n(m_rnd, n, dst);
                m_prnd = m_rnd + n;
            }
            else {
                m_prnd = endRnd();
            }
        }
        else {
            m_prnd = endRnd();
        }
    }

public:

    // constructors
    VMT19937Base()
        : m_prnd(nullptr)
        , m_prndEnd(m_rnd + s_n32InRndCache)
        , m_pst(nullptr)
        , m_pstEnd(m_state + s_N * s_n32inReg)
    {}
};

} // namespace Details
