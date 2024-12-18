#ifndef _RANDOM_H
#define _RANDOM_H

#include <inttypes.h>
#include <math.h>
#include "util.h"
#include "config.h"

#ifdef HAVE_SSE2
#include "SFMT-src-1.4/SFMT.h"
#else
#include "mt19937-64.h"
#endif

/*
 * Random number generator class.  Call seed() to seed, and next_uint64() or
 * next_double() to get random numbers.  search_location() will return a random
 * index from a cumulative probability array; rnorm() will return a number from
 * a normal distribution.
 *
 * The default constructor does not initialize the random state.  You MUST call
 * seed() before using the random number generator.
 *
 * This class contains the random state, so multiple threads can use this class,
 * provided that each one has its own instance seeded by a unique seed.
 *
 * There are 2 implementations for actual random number generation:  The normal
 * Mersenne twister, and the SSE2-accelerated Mersenne twister.  See the code
 * for details on those algorithms.
 */
class Random {
private:
#ifdef HAVE_SSE2
	sfmt_t state;
#else
	uint64_t state[MERSENNE_TWISTER_STATE_LEN];
#endif

public:
	Random() {
	}

	void seed(uint64_t seed) {
	#ifdef HAVE_SSE2
		sfmt_init_gen_rand(&state, seed);
	#else
		init_genrand64(seed, state);
	#endif
	}

	Random(uint64_t _seed) {
		seed(_seed);
	}
	uint64_t next_uint64() {
	#ifdef HAVE_SSE2
		return sfmt_genrand_uint64(&state);
	#else
		return genrand64_int64(state);
	#endif
	}
	double next_double()
	{
	#ifdef HAVE_SSE2
		return sfmt_genrand_res53(&state);
	#else
		return genrand64_real1(state);
	#endif
	}

	char next_base()
	{
		static const char bases[] = {'A', 'T', 'G', 'C'};
		return bases[next_uint64() % 4];
	}

	void random_dna_seq(char *seq, size_t len)
	{
		while (len--)
			*seq++ = next_base();
	}

	double rnorm(double mean, double sd)
	{
		double x1, x2, radius_pow2, y1 = 0.0;
		do {
			x1 = 2.0 * next_double() - 1.0;
			x2 = 2.0 * next_double() - 1.0;
			radius_pow2 = x1 * x1 + x2 * x2;
		} while (radius_pow2 >= 1.0 || radius_pow2 == 0);
		
		radius_pow2 = sqrt(-2.0 * log(radius_pow2) / radius_pow2);
		y1 = x1 * radius_pow2;

		return mean + y1 * sd;
	}

	// 
	// Binary search a 64-bit unsigned integer cumulative probability
	// distribution array.
	//
	// We are looking for the first location in @array of length @len that
	// has a value less than or equal to the selected random number.  The
	// index of this location is returned.
	//
	// @len must be nonzero, and @array must be a cumulative probability
	// distribution ending in the maximum value (2**64 - 1).
	//
	// Inline because this function is time critical (called at least once
	// for each called base pair) and it benefits from being inlined.
	size_t search_location(const uint64_t array[], size_t len) {
		uint64_t r = next_uint64();
		size_t left = 0;
		size_t right = len - 1;
		size_t middle;
		while (1) {
			middle = (left + right) >> 1;
			if (r > array[middle]) {
				left = middle + 1;
			} else {
				if (middle == 0 || r > array[middle - 1])
					return middle;
				right = middle - 1;
			}
		}
	}
};

class RandomBitGenerator : public Random {
private:
	unsigned bits_remaining;
	uint64_t bits;
public:
	RandomBitGenerator(uint64_t seed)
		: Random(seed), bits_remaining(0), bits(0)
	{ }
	unsigned next_bit() {
		unsigned ret;
		if (bits_remaining == 0) {
			bits = next_uint64();
			bits_remaining = 64;
		}
		bits_remaining--;
		ret = (bits & 1);
		bits >>= 1;
		return ret;
	}
};

#endif /* _RANDOM_H */
