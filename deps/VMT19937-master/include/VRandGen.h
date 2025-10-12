#pragma once

#include "VMT19937.h"
#include "VSFMT19937.h"

enum VRandGenQueryMode { QM_Any, QM_Scalar, QM_Block16, QM_StateSize };

namespace Details
{

template <typename GenBase, VRandGenQueryMode QueryMode>
class VRandGen : public GenBase
{
    typedef GenBase base_t;
public:
    typedef typename base_t::matrix_t matrix_t;
private:
    void completeStateInitialization(size_t nCommonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        // temporary workspace matrix
        BinaryMatrix<2, base_t::s_nMatrixBits> tmp;

        // apply common jump to state-0
        if (nCommonJumpRepeat) {
            MYASSERT(commonJump, "commonJump is required when nCommonJumpRepeat>0");

            // extract state 0
            base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

            for (size_t i = 0; i < nCommonJumpRepeat; ++i)
                commonJump->multiplyByColumn(tmp.rowBegin((i + 1) % 2), tmp.rowBegin(i % 2));

            // copy to the state vector shifting all bits to the right by 31
            base_t::vectorToState(0, (const uint32_t*)tmp.rowBegin(nCommonJumpRepeat % 2));
        }

        // if there are multiple states, distance them using the sequentialJump matrix
        if constexpr (base_t::s_nStates > 1) {
            if (sequentialJump) {
                // perform jump ahead of the s_regLenWords states
                // State_0 = State_0
                // State_1 = Jump x State_0
                // State_2 = Jump x State_1
                // ...

                // copy state 0 to the first row
                base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

                for (size_t s = 1; s < base_t::s_nStates; ++s) {
                    // multiply all rows by state s and store the result in pres

                    const uint8_t* psrc = (uint8_t*)tmp.rowBegin((s + 1) % 2);
                    uint8_t* pdst = (uint8_t*)tmp.rowBegin(s % 2);

                    sequentialJump->multiplyByColumn(pdst, psrc);

                    // copy to the state vector s
                    base_t::vectorToState(s, (const uint32_t*)pdst);
                }
            }
            else {
#if (VRANDGEN_TESTING!=1)
                THROW("Having multiple states and no sequential jump matrix does not make sense");
#endif
                // Copy state 0 to all other states.
                // Effectively the generator will produce s_nStates copies of each number
                // This only make sense for debugging purpose
                for (size_t w = 0; w < base_t::s_N; ++w)
                    for (size_t j = 0; j < base_t::s_n32InOneWord; ++j)
                        for (size_t s = 1; s < base_t::s_nStates; ++s)
                            base_t::m_state[w * base_t::s_n32inReg + s * base_t::s_n32InOneWord + j] = base_t::m_state[w * base_t::s_n32inReg + j];
            }
        }
    }


public:
    static const VRandGenQueryMode s_queryMode = QueryMode;

    VRandGen() {}

    VRandGen(uint32_t seed, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    VRandGen(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seeds, n_seeds, commonJumpRepeat, commonJump, sequentialJump);
    }

    // initializes m_state[s_N] with a seed
    void reinit(uint32_t s, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(s);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    void reinit(const uint32_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(seeds, nSeeds);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    // generates a random number on [0,0xffffffff] interval
    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(QueryMode == QM_Scalar || QueryMode == QM_Any);
        return base_t::genrand_uint32();
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QueryMode == QM_Block16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_uint32_blk16(dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_uint32_stateBlk(uint32_t* dst)
    {
        static_assert(QueryMode == QM_StateSize, "This function can only be invoked when query mode is QM_StateSize");
        base_t::genrand_uint32_stateBlk(dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(QueryMode == QM_Any, "This function can only be invoked when query mode is QM_Any");
        base_t::genrand_uint32_anySize(dst, n);
    }
};

} // namespace Details

template < size_t RegisterBitLen = SIMD_N_BITS
         , VRandGenQueryMode QueryMode = QM_Any
         , size_t RegisterBitLenImpl = std::min<size_t>(SIMD_N_BITS, RegisterBitLen)
         >
using VMT19937 = Details::VRandGen<Details::VMT19937Base<RegisterBitLen, RegisterBitLenImpl>, QueryMode>;

template < size_t RegisterBitLen = SIMD_N_BITS
         , VRandGenQueryMode QueryMode = QM_Any
         , size_t RegisterBitLenImpl = std::min<size_t>(SIMD_N_BITS, RegisterBitLen)
         >
using VSFMT19937 = Details::VRandGen<Details::VSFMT19937Base<RegisterBitLen, RegisterBitLenImpl>, QueryMode>;
