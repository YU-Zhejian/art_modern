#ifndef _READ_H
#define _READ_H

#include <vector>
#include <stddef.h>
#include <inttypes.h>
#include "config.h"

using std::vector;

class ReadPair;

/* Structure to hold information about a read indel.
 *
 * @ref_idx is the 0-based index in the reference read that the indel occurs
 * at; len is the length of the indel (positive for insertions, negative for
 * deletions.
 */
class Indel {
public:
	int ref_idx;
	int len;
	Indel(int _ref_idx, int _len) : ref_idx(_ref_idx), len(_len) { }
};

/*
 * Structure to represent a paired-end read.
 *
 * @pair is a reference to the containing ReadPair.
 */
class Read {
public:
	vector<char>   seq;
	vector<char>   raw_read;
	vector<char>   ref_read;
	vector<char>   quality_vals;
	vector<Indel>  indels;
	vector<int>    error_pos;
	ReadPair      &pair;
	int            mask_end_len;

	Read(ReadPair &_pair)
		: pair(_pair), mask_end_len(0)
	{ }

	inline int num_in_pair() const;
	inline char orientation() const;
};

/*
 * Structure to represent a pair of paired-end reads.
 */
class ReadPair {
public:
	Read         read_1;
	Read         read_2;
	const char  *ref_seq_id;
	const char  *ref_filename;
	int          insert_len;
	size_t       ref_seq_pos;
	uint64_t     pair_number;
	int          insert_len_mean;
	int          quality_shift;
	bool         reverse_order;
	bool         cyclicized;
	bool         dump_in_name;
	string       indiv_name;

	ReadPair()
		: read_1(*this), read_2(*this)
	{ }
};

inline int Read::num_in_pair() const
{
	if (this == &pair.read_1)
		return (pair.reverse_order) ? 2 : 1;
	else
		return (pair.reverse_order) ? 1 : 2;
}

inline char Read::orientation() const
{
	return (pair.reverse_order ^ pair.cyclicized ^
			(num_in_pair() == 1)) ? '+' : '-';
}


/* In the multi-threaded implementation, we pass Reads and ReadPairs between
 * threads in chunks.  Each chunk is a ReadSet or ReadPairSet, and each one
 * contains READS_PER_SET Reads or ReadPairs, respectively.
 *
 * The ReadSet merely contains pointers to the Reads, since these are changed
 * around, but the ReadPairSet contains the ReadPairs themselves.
 */
#ifdef ENABLE_THREADS
#define READS_PER_SET 16
#include "Lock.h"
class ReadPairSet {
private:
	unsigned nrefs;
	Lock nrefs_lock;
public:
	ReadPair *pairs;
	ReadPairSet() {
		pairs = new ReadPair[READS_PER_SET];
	}
	~ReadPairSet() {
		delete[] pairs;
	}

	// Set the reference count of the ReadPairSet
	void set_refs(unsigned n) {
		nrefs = n;
	}

	// Decrement the reference count of the ReadPairSet, returning true if
	// it's now 0.
	bool put_ref() {
		nrefs_lock.lock();
		unsigned n = --nrefs;
		nrefs_lock.unlock();
		return n == 0;
	}
};

class ReadSet {
public:
	Read **reads;
	ReadPairSet *pair_set;
	bool is_last;

	ReadSet() {
		reads = new Read*[READS_PER_SET];
	}
	~ReadSet() {
		delete[] reads;
	}
};
#endif // ENABLE_THREADS

#endif // _READ_H
