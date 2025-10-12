#ifndef GJRAND_HPP_INCLUDED
#define GJRAND_HPP_INCLUDED 1

/*
 * A C++ implementation of David Blackman's GJrand PRNG(s)
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018 Melissa E. O'Neill
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

/* Based on code provided by David Blackman, adapted for C++ */

#include <cstdint>

namespace gjrand_detail {

template <typename itype, typename rtype,
          itype c_init, itype d_inc,
          unsigned int p, unsigned int q, unsigned int r>
class gjrand {
protected:
    itype a_, b_, c_, d_;

    static constexpr unsigned int ITYPE_BITS = 8*sizeof(itype);
    static constexpr unsigned int RTYPE_BITS = 8*sizeof(rtype);

    static itype rotate(itype x, unsigned int k)
    {
        return (x << k) | (x >> (ITYPE_BITS - k));
    }

public:
    using result_type = rtype;
    using state_type = itype;

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return ~ result_type(0); }

    gjrand(itype seed1 = itype(0xcafef00dbeef5eedULL), itype seed2 = itype(0))
        : a_(seed1), b_(seed2), c_(c_init), d_(itype(0))
    {
        for (unsigned int j=14; j != 0; --j)
            advance();
    }

    void advance()
    {
        b_ += c_;
        a_ =  rotate(a_, p);
        c_ ^= b_;
        d_ += d_inc;
        a_ += b_;
        c_ =  rotate(c_, q);
        b_ ^= a_;
        a_ += c_;
        b_ =  rotate(b_, r);
        c_ += a_;
        b_ += d_;
    }
        

    rtype operator()()
    {
        advance();
        return rtype(a_);
    }

    bool operator==(const gjrand& rhs)
    {
        return (a_ == rhs.a_) && (b_ == rhs.b_) 
            && (c_ == rhs.c_) && (d_ == rhs.d_);
    }

    bool operator!=(const gjrand& rhs)
    {
        return !operator==(rhs);
    }

    // Not (yet) implemented:
    //   - arbitrary jumpahead (doable, but annoying to write).
    //   - I/O
    //   - Seeding from a seed_seq.
};

} // end namespace

///// ---- Specific GJrand Generators ---- ////
//
// Each size has variations corresponding to different parameter sets.
// Each variant will create a distinct (and hopefully statistically
// independent) sequence.
//

// - 256 state bits, uint64_t output
//   This version has been extensively tested by David Blackman.

using gjrand64 =
    gjrand_detail::gjrand<uint64_t, uint64_t, 2000001, 0x55aa96a5, 32,23,19>;

// - 128 state bits, uint32_t output
//   These parameters are from David Blackman, but he has *NOT*
//   extensively tested it.

using gjrand32 =
    gjrand_detail::gjrand<uint32_t, uint32_t, 2000001, 0x96a5, 16, 11, 19>;

// TINY VERSIONS FOR TESTING AND SPECIALIZED USES ONLY

// - 64 state bits, uint16_t output

using gjrand16 =
    gjrand_detail::gjrand<uint16_t, uint16_t, 2001, 0x96a5, 8, 5, 10>;

// - 32 state bits, uint8_t output

using gjrand8 =
    gjrand_detail::gjrand<uint8_t, uint8_t, 201, 0x35, 4, 2, 5>;

#endif // GJRAND_HPP_INCLUDED
