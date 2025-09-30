#pragma once

#include "SIMD.h"
#include "codecs.h"
#include "utils.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>

inline constexpr uint8_t bitmask(size_t b)
{
    return (uint8_t)(uint8_t(1) << b);
}

// simd implementation
template <size_t SimdBits>
struct BinaryVectorMultiplier
{
    // multyply the row pointed by pr by all the nRows pointed by pc
    // all rows must have size nBitCols and be padded to the next multiple of the simd vector length
    template <size_t nRows, size_t nBitCols, size_t nBitColsPadded>
    static uint8_t multiply8(const uint8_t* _pr, const uint8_t** _pc)
    {
        typedef Details::SimdRegister<SimdBits, SimdBits> simd_t;

        const size_t nSimdBytes = sizeof(simd_t);
        const size_t nSimdBits = 8 * nSimdBytes;
        const size_t nSimdStep = nBitCols / nSimdBits + ((nBitCols % nSimdBits) != 0);
        static_assert(nSimdStep * nSimdBits <= nBitColsPadded);

        simd_t accumulator[nRows];
        for (auto& a : accumulator)
            a = simd_t::zero();
        for (size_t i = 0; i < nSimdStep; ++i) {
            simd_t a(_pr + i * nSimdBytes);
            for (size_t c = 0; c < nRows; ++c) {
                simd_t b(_pc[c] + i * nSimdBytes);
                simd_t p = a & b;
                accumulator[c] = accumulator[c] ^ p;
            }
        }

        uint8_t result = 0;
        static_assert(nRows <= 8 * sizeof(result));

        for (size_t c = 0; c < nRows; ++c)
            result |= accumulator[c].parity() << c;

        return result;
    }
};


template <size_t _nBitRows, size_t _nBitCols = _nBitRows>
class BinaryMatrix
{
    static const size_t s_nAlignBits = 512;  // suitable for AVX512
    static const size_t s_nAlignBytes = s_nAlignBits / 8;  // suitable for AVX512
    static_assert(s_nAlignBits == 8 * s_nAlignBytes);

    typedef BinaryMatrix<_nBitRows, _nBitCols> this_t;

    AlignedVector<uint8_t, s_nAlignBytes> m_data;

public:
    static const size_t s_nBitRows = _nBitRows;
    static const size_t s_nBitCols = _nBitCols;
    static const size_t s_nPaddingBits = (s_nAlignBits - (s_nBitCols % s_nAlignBits)) % s_nAlignBits;
    static const size_t s_nBitColsPadded = s_nBitCols + s_nPaddingBits;
    static const size_t s_nBytesPerRow = s_nBitCols / 8+ (s_nBitCols % 8 != 0);
    static const size_t s_nBytesPerPaddedRow = s_nBitColsPadded / 8;  // inclusive of padding
    static const size_t s_nUsedBytes = s_nBitRows * s_nBytesPerPaddedRow;
    static const size_t s_binStreamSize = s_nBitRows * s_nBytesPerRow;

private:

    struct OMatrixStream : public std::stringbuf
    {
        this_t* m;
        size_t rowIndex;
        char* ptr;
        const char* colend;
        OMatrixStream(this_t* matrix) : m(matrix), rowIndex(size_t(-1)), ptr(nullptr), colend(nullptr) {}
        virtual std::stringbuf::int_type overflow(std::stringbuf::int_type ch)
        {
            if (ptr != colend) {
                *ptr++ = (char)ch;
                return 0;
            }
            if (++rowIndex < s_nBitRows) {
                ptr = (char*)m->rowBegin(rowIndex);
                colend = ptr + s_nBytesPerRow;
                *ptr++ = (char)ch;
                return 0;
            }
            return std::char_traits<char>::eof();
        }
    };


    uint8_t& getByte(size_t rowIndex, size_t colByteIndex) { return m_data[rowIndex * s_nBytesPerPaddedRow + colByteIndex]; }
    const uint8_t& getByte(size_t rowIndex, size_t colByteIndex) const { return m_data[rowIndex * s_nBytesPerPaddedRow + colByteIndex]; }

    // compute the index in the uint8_t matrix and the offset within the uint8_t element
    static auto bitPos(size_t bitColIndex)
    {
        size_t byteIndex = bitColIndex / 8;
        uint8_t bitOffset = bitColIndex % 8;
        return std::make_pair(byteIndex, bitOffset);
    }

public:
    BinaryMatrix(bool allocateMemory = true)
    {
        if (allocateMemory) {
            // allocate memory
            m_data.init(s_nUsedBytes);
            // initialize to zero
            resetZero();
        }
    }

    bool operator==(const BinaryMatrix& rhs) const
    {
        if (0 == memcmp(m_data.data(), rhs.m_data.data(), s_nUsedBytes))
            return true;
#if (VRANDGEN_TESTING==1)
        for (size_t i = 0; i < s_nBitRows; ++i) {
            for (size_t j = 0; j < s_nBitCols; ++j) {
                if (getBit(i, j) != rhs.getBit(i, j)) {
                    std::cout << "diff at: (" << i << "," << j << "): " << getBit(i, j) << " vs " << rhs.getBit(i, j) << "\n";
                    break;
                }
            }
        }
#endif
        return false;
    }

    // set bit at row r and column c
    // r and c are bit indices
    void setBit(size_t bitRowIndex, size_t bitColIndex)
    {
        auto [byteIndex, bitOffset] = bitPos(bitColIndex);
        getByte(bitRowIndex, byteIndex) |= bitmask(bitOffset);
    }

    // get bit at row r and column c
    // r and c are bit indices
    bool getBit(size_t bitRowIndex, size_t bitColIndex) const
    {
        auto [byteIndex, bitOffset] = bitPos(bitColIndex);
        return getByte(bitRowIndex, byteIndex) & bitmask(bitOffset);
    }

    void printBits(size_t bitRowIndex, size_t bitColIndex, size_t nRows, size_t nCols, bool matlabIndexStyle = true) const
    {
        std::cout << "bit columns " << bitColIndex + matlabIndexStyle << " to " << bitColIndex + nCols - (1-matlabIndexStyle) << "\n";
        for (size_t i = 0; i < nRows; ++i) {
            std::cout << std::setw(5) << bitRowIndex + i + matlabIndexStyle << ": ";
            for (size_t j = 0; j < nCols; ++j)
                std::cout << getBit(bitRowIndex + i, bitColIndex + j);
            std::cout << std::endl;
        }
    }

    // count non-zero bit in a range
    size_t nnz(size_t bitRowIndex, size_t bitColIndex, size_t nRows, size_t nCols) const
    {
        size_t n = 0;
        for (size_t i = 0; i < nRows; ++i) {
            for (size_t j = 0; j < nCols; ++j)
                n += getBit(bitRowIndex + i, bitColIndex + j);
        }
        return n;
    }

    void toMatlab(const char *path, const char *name) const
    {
        std::ostringstream fn;
        fn << path << name << ".m";
        std::ofstream os(fn.str());
        std::vector<size_t> I, J, Z;
        os << "function y=" << name << "()\n";
        for (size_t i = 0; i < s_nBitRows; ++i) {
            for (size_t j = 0; j < s_nBitCols; ++j)
                if (getBit(i, j)) {
                    I.push_back(i + 1);
                    J.push_back(j + 1);
                    Z.push_back(1);
                }
        }
        os << "I=[";
        for (auto i : I) os << i << " ";
        os << "];\n";
        os << "J=[";
        for (auto j : J) os << j << " ";
        os << "];\n";
        os << "Z=[";
        for (auto z : Z) os << z << "==1 ";
        os << "];\n";
        os << "y = sparse(I, J, Z, " << s_nBitRows << ", " << s_nBitCols << ");\n";
        os << "end\n";
    }

    template <typename OS>
    void txtRowEncoder(OS& os, void (*encoder)(std::string&, const std::string&)) const
    {
        // FIXME: too many copies
        std::string binStr;
        binStr.resize(s_nBytesPerRow * s_nBitRows);
        char* pout = &binStr[0];
        for (size_t r = 0; r < s_nBitRows; ++r) {
            auto prow = (const char*)rowBegin(r);
            std::copy(prow, prow + s_nBytesPerRow, pout);
            pout += s_nBytesPerRow;
        }
        
        std::string encodedStr;
        (*encoder)(encodedStr, binStr);
        os << encodedStr;
    }

    template <typename IS>
    void txtRowDecoder(IS& is, void (*decoder)(std::string&, const std::string&))
    {
        std::string encodedStr, binStr;
        std::getline(is, encodedStr);
        (*decoder)(binStr, encodedStr);

        MYASSERT(binStr.length() == s_nBytesPerRow * s_nBitRows, "bad stream length, expect " << (s_nBytesPerRow * s_nBitRows) << ", got " << binStr.length() << "\n");

        const uint8_t* p = (const uint8_t*)&binStr[0];
        for (size_t r = 0; r < s_nBitRows; ++r, p += s_nBytesPerRow)
            std::copy(p, p + s_nBytesPerRow, (char *) rowBegin(r));
    }

    void txtRowDecoderStream(std::istream& is, std::ostream& (*decoder)(std::ostream&, std::istream&))
    {
        OMatrixStream buffer(this);
        std::ostream os(&buffer);
        (*decoder)(os, is);
    }

    template <typename OS>
    void toBin(OS& os) const
    {
        for (size_t r = 0; r < s_nBitRows; ++r)
            os.write((const char *) rowBegin(r), s_nBytesPerRow);
    }

    template <typename IS>
    void fromBin(IS& is)
    {
        for (size_t r = 0; r < s_nBitRows; ++r)
            is.read((char*) rowBegin(r), s_nBytesPerRow);
    }


    void fromArrayChar(const uint8_t*pchar, size_t len)
    {
        MYASSERT(len == s_nBitRows * s_nBytesPerRow, "array length must be same length" << s_nBitRows * s_nBytesPerRow);
        for (size_t r = 0; r < s_nBitRows; ++r)
            std::copy_n(pchar + s_nBytesPerRow * r, s_nBitRows * s_nBytesPerRow, (uint8_t*)rowBegin(r));
    }

    template <typename OS>
    void toBase64(OS& os) const
    {
        txtRowEncoder(os, &Encoder::textToBase64);
    }

    template <typename IS>
    void fromBase64(IS& is)
    {
        //txtRowDecoder(os, &Encoder::base64ToText);
        txtRowDecoderStream(is, &Encoder::base64ToTextStream);
    }


    template <typename OS>
    void toHex(OS& os) const
    {
        txtRowEncoder(os, &Encoder::textToHex);
    }

    template <typename IS>
    void fromHex(IS& is)
    {
        //txtRowDecoder(is, &Encoder::hexToText);
        txtRowDecoderStream(is, &Encoder::hexToTextStream);
    }

    template <typename OS>
    void toArrayChar(OS& os) const
    {
        for (size_t r = 0; r < s_nBitRows; ++r) {
            const char* p = (const char*)rowBegin(r);
            for (size_t c = 0; c < s_nBytesPerRow; ++c)
                os << (unsigned) (uint8_t) p[c] << ',';
        }
    }

    size_t nnz() const
    {
        size_t n = 0;
        static_assert(s_nUsedBytes % sizeof(uint64_t) == 0);
        for (const uint64_t *p = (const uint64_t*)m_data.data(), * const pend = p + (s_nUsedBytes / sizeof(*p)); p != pend; ++p)
            n += popcnt(*p);
        return n;
    }

    auto sparseIndex() const
    {
        size_t n = nnz();
        return std::make_pair(n, double(n) / (double(s_nBitCols) * double(s_nBitRows)));
    }

    // returns the word of type U containing the bit in position colBitIndex
    template <typename U>
    U getWordU(size_t rowBitIndex, size_t colBitIndex) const
    {
        // requires that colBitIndex is a multiple of sizeof(T)
        auto p = (const U*)rowBegin(rowBitIndex);
        return p[colBitIndex / (8 * sizeof(U))];
    }

    const uint8_t* rowBegin(size_t rowIndex) const { return m_data.data() + rowIndex * s_nBytesPerPaddedRow; }
    uint8_t* rowBegin(size_t rowIndex) { return m_data.data() + rowIndex * s_nBytesPerPaddedRow; }

    void resetZero() { std::fill_n(m_data.data(), s_nUsedBytes, 0); }

    // Multiply all rows by 1 column (psrc) and stores the resulting column in (pdst)
    // It is assumed that:
    //   psrc and dst have at least size s_nBytesPerPaddedRow
    //   psrc contains zeros in the padded area
    //   pdst contains all zeros
    void multiplyByColumn(uint8_t* pdst, const uint8_t* psrc) const
    {
        const uint8_t* rowptrs[8];
        size_t r;
        for (r = 0; r + 7 < s_nBitRows; r += 8) {
            for (size_t h = 0; h < 8; ++h)
                rowptrs[h] = rowBegin(r + h);
            pdst[r / 8] = BinaryVectorMultiplier<SIMD_N_BITS>::template multiply8<8, s_nBitCols, s_nBitColsPadded>(psrc, rowptrs);
        }
        const size_t nResidualRows = s_nBitRows % 8;
        if constexpr (nResidualRows)
        { // this takes care of the residual rows, in case it is not a multiple of 8
            for (size_t h = 0; h < nResidualRows; ++h)
                rowptrs[h] = rowBegin(r + h);
            pdst[r / 8] = BinaryVectorMultiplier<SIMD_N_BITS>::template multiply8<nResidualRows, s_nBitCols, s_nBitColsPadded>(psrc, rowptrs);
        }
    }

    // initialize from a file saved in matlab with the commands:
    //    [I,J,V] = find(sparse(M));
    //    U = [I,J];
    //    save('filename','U', '-ascii');
    void fromMatlabSparseFile(const char* filename)
    {
        std::ifstream is(filename);
        if (!is.good()) {
            std::cout << "failure opening " << filename << "\n";
            std::exit(-1);
        }
        double r, c;
        while (true) {
            is >> r >> c;
            if (!is.eof())
                setBit((size_t)r - 1, (size_t)c - 1);
            else
                break;
        }
    }

    void initRand()
    {
        for (size_t r = 0; r < s_nBitRows; ++r) {
            //uint8_t* pr = rowBegin(r);
            size_t c = 0;
            {
                typedef uint16_t word_t;
                for (; c + sizeof(word_t) * 8 < s_nBitCols; c += sizeof(word_t) * 8) {
                    word_t rnd = word_t(rand()) & word_t(-1);
                    ((word_t*)rowBegin(r))[c / (sizeof(word_t) * 8)] = rnd;
                }
            }
            {
                typedef uint8_t word_t;
                for (; c + sizeof(word_t) * 8 < s_nBitCols; c += sizeof(word_t) * 8) {
                    word_t rnd = word_t(rand()) & word_t(-1);
                    ((word_t*)rowBegin(r))[c / (sizeof(word_t) * 8)] = rnd;
                }
            }
            for (;  c < s_nBitCols; ++c) {
                auto rnd = rand() % 2;
                if (rnd)
                    setBit(r, c);
            }
        }
    }

};

