/*
 * util.cpp
 *
 * Various useful functions.
 */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vector>
using std::vector;

#include "util.h"
#include "InputStream.h"

#define LOG_TAG "[pIRS] "

void fatal_error(const char *msg, ...)
{
	fputs("************************************************************\n", stderr);
	va_list va;
	va_start(va, msg);
	fflush(stdout);
	fputs(LOG_TAG "FATAL ERROR:\n", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	fputs("************************************************************\n", stderr);
	exit(1);
}

void fatal_error_with_errno(const char *msg, ...)
{
	fputs("************************************************************\n", stderr);
	va_list va;
	va_start(va, msg);
	int errno_save = errno;
	fflush(stdout);
	fputs(LOG_TAG "FATAL ERROR:\n", stderr);
	vfprintf(stderr, msg, va);
	fprintf(stderr, ": %s\n", strerror(errno_save));
	va_end(va);
	fputs("************************************************************\n", stderr);
	exit(1);
}

FILE *info_messages_fp = stdout;

void info(const char *format, va_list va)
{
	if (info_messages_fp) {
		fputs(LOG_TAG, info_messages_fp);
		vfprintf(info_messages_fp, format, va);
		fflush(info_messages_fp);
	}
}

void info(const char *format, ...)
{
	va_list va;
	va_start(va, format);
	info(format, va);
	va_end(va);
}

void warning(const char *format, ...)
{
	va_list va;
	va_start(va, format);
	fputs("************************************************************\n"
	      LOG_TAG"WARNING:\n", stderr);
	vfprintf(stderr, format, va);
	fputs("\n************************************************************\n", stderr);
	va_end(va);
}

//from ASCII of A C G T to 0 1 2 3, auto dealing with upper or lower case.
//8bit char type, A=a=0, C=c=1, G=g=2, T=t=3, others as 4.
const char dna_char_to_bin_tab[256] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,

	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

const char dna_bin_to_char_tab[5] = {
	'A', 'C', 'G', 'T', 'N',
};


const char dna_char_complement_tab[256] = {
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'

};

/* Reads the next scaffold from a FASTA file and provides the id and sequence as
 * strings.  Returns %true on success, %false on EOF, and exits on error. */
bool read_scaffold(InputStream &in, string &id, vector<char> &seq)
{
	bool ret = false;
	char *line = NULL;
	size_t n;
	ssize_t len;
	int c;


	if ((len = in.getline(&line, &n)) == -1)
		goto out;

	if (line[0] != '>')
		fatal_error("Expected FASTA tag line in reference sequence "
			    "\"%s\"", in.s_filename);
	trim(line + 1, len - 1);
	id = line + 1;

	info("Reading scaffold \"%s\" into memory\n", id.c_str());

	ret = true;
	seq.clear();
	while (1) {
		c = in.getc();
		if (c == -1)
			goto out;
		in.ungetc(c);
		if (c == '>')
			goto out;
		if ((len = in.getline(&line, &n)) == -1)
			goto out;
		len = trim(line, len);
		seq.insert(seq.end(), line, line + len);
	}
out:
	free(line);
	return ret;
}

char *xstrdup(const char *s)
{
	size_t len = strlen(s);
	char *p = (char*)xmalloc(len + 1);
	return (char*)memcpy(p, s, len + 1);
}

char *xstrdup2(const char *s1, const char *s2)
{
	size_t len1 = strlen(s1);
	size_t len2 = strlen(s2);
	size_t len = len1 + len2;
	char *p = (char*)xmalloc(len + 1);
	memcpy(p, s1, len1);
	memcpy(p + len1, s2, len2 + 1);
	return p;
}

void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		fatal_error("Out of memory: tried to allocate %lu bytes", size);
	return p;
}

void *xrealloc(void *ptr, size_t size)
{
	void *p = realloc(ptr, size);
	if (!p)
		fatal_error("Out of memory: tried to reallocate %lu bytes", size);
	return p;
}

static void *new_matrix_dimv(size_t dimc, const size_t dimv[])
{
	if (dimc > 1) {
		uint64_t **p = (uint64_t**)xmalloc(dimv[0] * sizeof(uint64_t*));
		for (size_t i = 0; i < dimv[0]; i++)
			p[i] = (uint64_t*)new_matrix_dimv(dimc - 1, dimv + 1);
		return p;
	} else {
		uint64_t *p = (uint64_t*)xmalloc(dimv[0] * sizeof(uint64_t));
		zero_array(p, dimv[0]);
		return p;
	}
}

static void delete_matrix_dimv(void *matrix, size_t dimc, const size_t dimv[])
{
	if (matrix) {
		if (dimc > 1)
			for (size_t i = 0; i < dimc; i++)
				delete_matrix_dimv(((void**)matrix)[i], 
						   dimc - 1, dimv + 1);
		free(matrix);
	}
}

// Allocate a @dimc-dimensional matrix of uint64_t's.  The variadic parameters
// provide the dimensions (there must be exactly @dimc such parameters, and they
// must all be greater than 0.)
void *new_matrix(size_t dimc, ...)
{
	va_list va;
	va_start(va, dimc);
	size_t dimv[dimc];
	for (size_t i = 0; i < dimc; i++)
		dimv[i] = va_arg(va, size_t);
	return new_matrix_dimv(dimc, dimv);
}

//
// Delete a multi-dimensional matrix of 64-bit unsigned integers.
//
void delete_matrix(void *matrix, size_t dimc, ...)
{
	va_list va;
	va_start(va, dimc);
	size_t dimv[dimc];
	for (size_t i = 0; i < dimc; i++)
		dimv[i] = va_arg(va, size_t);
	return delete_matrix_dimv(matrix, dimc, dimv);
}

//
// Prepares a cumulative probability array.
//
// The array is made such that we can generate a random 64-bit integer and
// search for the first array element that is less than or equal to the random
// number, and that array slot tells us something (such as indel length, quality
// score, base pair etc. depending on where the array is being used).
//
// 64-bit unsigned integers are used instead of doubles for better efficiency
// and precision.  A probability of 1 is 2**64 - 1, or 0xffffffffffffffff.
//
// The array is an input-output parameter for this function.  It initially
// contains a list of raw probabilities that are neither cumulative nor scaled
// to 2**64-1.  This function will make the probability array cumulative; then it
// will re-scale the probabilities such that 1 is 2**64-1.  There may be a
// remainder; this remainder is assigned to the element at @extra_idx.  This
// must be done, so that last array element can contain 2**64-1.
//
bool prepare_uint64_probability_array(uint64_t array[], size_t len,
				      size_t extra_idx)
{
	if (len == 0)
		return true;
	uint64_t sum = array[0];
	for (size_t i = 1; i < len; i++) {
		if (sum + array[i] < sum)
			return false;
		sum += array[i];
	}
	if (sum == 0) {
		for (size_t i = 0; i < extra_idx; i++)
			array[i] = 0;
		for (size_t i = extra_idx; i < len; i++)
			array[i] = ~0ULL;
	} else {

		uint64_t scale_factor = ~0ULL / sum;
		uint64_t extra = ~0ULL % sum;

		for (size_t i = 0; i < len; i++)
			array[i] *= scale_factor;
		array[extra_idx] += extra;
		uint64_t cum_sum = 0;
		for (size_t i = 0; i < len; i++) {
			cum_sum += array[i];
			array[i] = cum_sum;
		}
	}
	return true;
}


//
// Reverse-complements a DNA sequence in place.
//
void reverse_complement(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p <= pp) {
		char tmp = *p;
		*p++ = dna_char_complement(*pp);
		*pp-- = dna_char_complement(tmp);
	}
}

//
// Trims whitespace from the end of a string and returns the new length.
size_t trim(char *s, size_t len)
{
	while (len != 0 && isspace(s[len - 1]))
		s[--len] = '\0';
	return len;
}
