#pragma once

#include "bit_matrix.h"
#include "Params.h"

#include <thread>
#include <atomic>
#include <memory>
#include <sstream>

template <size_t N>
struct BinarySquareMatrix : BinaryMatrix<N, N>
{
    typedef BinaryMatrix<N,N> base_t;

    typedef uint32_t word_t;

    static const size_t s_nBits = N;
    static const size_t s_nWordBits = sizeof(word_t) * 8;

    static const size_t s_nBitCols = base_t::s_nBitCols;
    static const size_t s_nBitColsPerBlk = 8;
    const size_t s_nColBlks = s_nBits / s_nBitColsPerBlk;
    typedef BinaryMatrix<s_nBitColsPerBlk, base_t::s_nBitRows> buffer_t;

    template <size_t NRows, size_t nColumns, size_t...Is>
    static void transpose16x8(buffer_t& buffer, const base_t& src, size_t rowBitIndex, size_t colBitIndex,  std::index_sequence<Is...>&&)
    {
        //const size_t simdBits = 128;
        //const size_t nRowsPerBlk = simdBits / s_nBitColsPerBlk;

        static_assert(s_nBitColsPerBlk == 8);

        const size_t n = sizeof...(Is) - 1;
        auto cs = _mm_set_epi8((n - Is < NRows ? src.template getWordU<uint8_t>(rowBitIndex + (n - Is), colBitIndex) : 0)...);
        
        for (size_t cb = 0; cb < nColumns; ++cb) {
            auto mask = _mm_movemask_epi8(_mm_slli_epi64(cs, (int)(s_nBitColsPerBlk - 1 - cb)));
            uint16_t* pc = (uint16_t*)buffer.rowBegin(cb);
            pc[rowBitIndex / (8 * sizeof(uint16_t))] = (uint16_t)mask;
        }
    }

    // transpose columns from colBitIndex to colBitIndex+nBitColsPerBlk-1 and copies then to colBuffer
    // so that then is it more efficient to compute the square of the matrix, where we have to multiply rows by columns
    template <size_t nColumns>
    static void transposeColumns(buffer_t& colBuffer, const base_t& src, size_t colBitIndex)
    {
        const size_t simdBits = 128;
        const size_t nRowsPerBlk = simdBits / s_nBitColsPerBlk;
        const size_t nRowBlks = s_nBits / nRowsPerBlk;
        const size_t nResidualRows = base_t::s_nBitRows - nRowBlks * nRowsPerBlk;

        colBuffer.resetZero();

        for (size_t i = 0; i < nRowBlks; ++i)
            transpose16x8<nRowsPerBlk, nColumns>(colBuffer, src, i * nRowsPerBlk, colBitIndex, std::make_index_sequence<nRowsPerBlk>{});

        if (nResidualRows)
            transpose16x8<nResidualRows, nColumns>(colBuffer, src, nRowBlks * nRowsPerBlk, colBitIndex, std::make_index_sequence<nRowsPerBlk>{});
    }

    template <size_t nColumns, typename DST>
    void multiplyBlock(DST& colBuffer, const base_t& src, size_t colBitIndex)
    {
        // transpose columns from colBitIndex to colBitIndex+nBitColsPerBlk-1 and copies then to colBuffer
        // so that then is it more efficient to compute the square of the matrix, where we have to multiply rows by columns
        transposeColumns<nColumns>(colBuffer, src, colBitIndex);

#if 0 || defined(_DEBUG)
        static bool firstError = true;
        for (size_t r = 0; r < base_t::s_nBitRows; ++r) {
            if (!firstError)
                break;
            for (size_t c = 0; c < nColumns; ++c)
                if (src.getBit(r, c + colBitIndex) != colBuffer.getBit(c, r)) {
                    std::cout << "error at coord " << r << "," << c << "\n";
                    src.printBits(r / 16 * 16, (c + colBitIndex) / 8 * 8, 16, 8);
                    colBuffer.printBits(0, r / 16 * 16, 8, 16);
                    firstError = true;
                    break;
                }
        }
#endif

        const uint8_t* pcs[nColumns];
        for (size_t c = 0; c < nColumns; ++c)
            pcs[c] = colBuffer.rowBegin(c);

        for (size_t r = 0; r < s_nBits; ++r) {
            auto pr = src.rowBegin(r);
            uint8_t active = BinaryVectorMultiplier<SIMD_N_BITS>::multiply8<nColumns, base_t::s_nBitCols, base_t::s_nBitColsPadded>(pr, pcs);
            this->rowBegin(r)[colBitIndex / 8] = active;
        }
    }

    void squareBlock(const base_t& src, buffer_t& colBuffer, std::atomic<size_t>& nextJobIndex)
    {
        while (true) {
            const size_t curJobIndex = nextJobIndex++;
            const size_t colBitIndex = curJobIndex * s_nBitColsPerBlk;
            if (curJobIndex < s_nColBlks) {
                multiplyBlock<s_nBitColsPerBlk>(colBuffer, src, colBitIndex);
            }
            else if (curJobIndex == s_nColBlks) { // in the last block we have less than s_nBitColsPerBlk columns to process
                constexpr size_t nColumns = s_nBitCols % s_nBitColsPerBlk;
                if constexpr (nColumns)
                    multiplyBlock<nColumns>(colBuffer, src, colBitIndex);
            }
            else
                break;
        }
    }

    void square(const base_t& src, std::vector<buffer_t>& buffers)
    {
        std::vector< std::shared_ptr<std::thread>> threads(buffers.size() - 1);

        std::atomic<size_t> nextJobIndex(0);

        // launch participating threads
        for (size_t i = 0; i < threads.size(); ++i)
            threads[i].reset(new std::thread(&BinarySquareMatrix::squareBlock, this, std::cref(src), std::ref(buffers[i]), std::ref(nextJobIndex)));

        // main thread participate to the loop
        squareBlock(src, buffers.back(), nextJobIndex);

        for (auto& th : threads)
            th->join();
    }

    void printSparsity()
    {
        auto [n, frac] = this->sparseIndex();
        std::cout << "nnz: " << n << ", sparse index: " << std::fixed << std::setprecision(4) << 100.0 * frac << "%\n";
    }
};

struct MT19937Matrix : BinarySquareMatrix<Details::MT19937Params::s_nMatrixBits>
{
    typedef BinarySquareMatrix<Details::MT19937Params::s_nMatrixBits> base_t;

    // Initialize the matrix as per MT19937 32 bit generator transition matrix
    // This is equivalent to a jump ahead of 1 random number
    MT19937Matrix()
    {
        init1();
    }

    MT19937Matrix(const std::string& binaryfilename)
    {
        fromBinaryFile(binaryfilename);
    }

    MT19937Matrix(const uint8_t* pchar, size_t len)
    {
        fromArrayChar(pchar, len);
    }

    template <size_t N>
    MT19937Matrix(const uint8_t(&pchar)[N])
    {
        fromArrayChar(pchar, N);
    }

    void init1()
    {
        static const size_t s_nBits = base_t::s_nBitRows;
        static const size_t s_nWordBits = base_t::s_nWordBits;
        static const uint32_t s_matA = Details::MT19937Params::s_matrixA;
        static const uint32_t s_M = Details::MT19937Params::s_M;

        // from row 0 to to row nBits - 32, state bits are just shifted left by 32 bits
        for (uint32_t r = 0; r < s_nBits - s_nWordBits; ++r)
            setBit(r, r + s_nWordBits);

        // the new state element is composed of element which was in position M
        for (uint32_t i = 0; i < s_nWordBits; i++)
            setBit(s_nBits - s_nWordBits + i, 1 + (s_M - 1) * s_nWordBits + i);
        for (uint32_t i = 0; i < s_nWordBits; ++i)
            if (s_matA & (uint32_t(1) << i))
                setBit(s_nBits - s_nWordBits + i, 1);
        setBit(s_nBits - 2, 0);
        for (uint32_t i = 0; i < s_nWordBits - 2; ++i)
            setBit(s_nBits - s_nWordBits + i, 2 + i);
    }

    // initialize from a binary file saved with the toBin method
    void fromBinaryFile(const std::string& filename)
    {
        std::ifstream is(filename, std::ios::binary);
        MYASSERT(is.is_open(), "error opening binary file: " << filename);
        base_t::fromBin(is);
#if (VRANDGEN_TESTING==1)
        std::cout << "loaded matrix from file: " << filename << "\n";
        printSparsity();
#endif
    }

    // initialize from a binary file saved with the toBin method
    void fromArrayChar(const uint8_t* pchar, size_t len)
    {
        base_t::fromArrayChar(pchar, len);
    }
};

struct SFMT19937Matrix : BinarySquareMatrix<Details::SFMT19937Params::s_nMatrixBits>
{
    typedef BinarySquareMatrix<Details::SFMT19937Params::s_nMatrixBits> base_t;

    // Initialize the matrix as per MT19937 32 bit generator transition matrix
    // This is equivalent to a jump ahead of 4 random numbers
    SFMT19937Matrix()
    {
        init4();
    }

    SFMT19937Matrix(const std::string& binaryfilename)
    {
        fromBinaryFile(binaryfilename);
    }

    void init4()
    {
    	using namespace Details;

        const uint32_t masks[] =
            { SFMT19937Params::s_SFMT_MSK1
            , SFMT19937Params::s_SFMT_MSK2
            , SFMT19937Params::s_SFMT_MSK3
            , SFMT19937Params::s_SFMT_MSK4
            };

        auto K = [&masks](size_t i) -> bool { return masks[i / 32] & (1 << (i % 32)); };

        const size_t s_nBits = base_t::s_nBitRows;
        const int s_N = SFMT19937Params::s_N;
        const int s_M = SFMT19937Params::s_M;

        // from row 0 to to row nBits - 128, state bits are just shifted left by 128 bits
        for (uint32_t r = 0; r < s_nBits - 128; ++r)
            setBit(r, r + 128);

        // The new state element is composed of elements which were in position {0, M, N-2, N-1}
        // W[N] = w[0] + (w[0] <<< 8) + ((w[M] >> 11) & K) + (w[N-2] >>> 8) + (w[N - 1] << 18)
        // where:
        //    <<< and >>> are 128 bit shifts
        //    << and >> are 32 bit shifts

        const size_t B = s_nBits - 128;

        // w[0]
        for (size_t r = 0; r < 128; ++r)
            setBit(B + r, r);

        // w[0] <<< 8
        for (size_t r = 8; r < 128; ++r)
            setBit(B + r, r - 8);

        // w[N-2] >>> 8
        for (size_t r = 0; r < 128 - 8; ++r)
            setBit(B + r, (s_N-2)*128 + r + 8);

        // w[N - 1] << 18
        for (size_t w = 0; w < 4; ++w)
            for (size_t r = 18; r < 32; ++r)
                setBit(B + 32 * w + r, (s_N - 1) * 128 + 32 * w + r - 18);

        // K & (w[M] >> 11)
        for (size_t w = 0; w < 4; ++w)
            for (size_t r = 0; r < 32-11; ++r)
                if (K(32 * w + r))
                    setBit(B + 32 * w + r, s_M * 128 + 32 * w + r + 11);
    }

    // initialize from a binary file saved with the toBin method
    void fromBinaryFile(const std::string& filename)
    {
        std::ifstream is(filename, std::ios::binary);
        MYASSERT(is.is_open(), "error opening binary file: " << filename);
        base_t::fromBin(is);
#if (VRANDGEN_TESTING==1)
        std::cout << "loaded matrix from file: " << filename << "\n";
        printSparsity();
#endif
    }
};
