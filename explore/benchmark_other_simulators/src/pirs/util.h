/*
 * util.h
 *
 * Various useful functions.
 */
#ifndef _UTIL_H
#define _UTIL_H

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
# 		define __cold __attribute__((cold))
#	else
#		define __cold
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#else
#	define __noreturn
#	define __cold
#	define __format
#endif

#include <stdio.h>
#include <inttypes.h>

#define GCC_VERSION (__GNUC__ * 10000 \
		+ __GNUC_MINOR__ * 100 \
		+ __GNUC_PATCHLEVEL__)
#if GCC_VERSION >= 40201
#	ifndef WITH_DEBUG
#		pragma GCC diagnostic ignored "-Wfloat-equal"
#	endif
#endif

class InputStream;

extern FILE *info_messages_fp;

extern void fatal_error(const char *msg, ...)
	       		__cold __noreturn __format(printf, 1, 2);

extern void fatal_error_with_errno(const char *msg, ...)
	       		__cold __noreturn __format(printf, 1, 2);

extern void info(const char *msg, ...)
	       		__cold __format(printf, 1, 2);

extern void info(const char *msg, va_list va) __cold __format(printf, 1, 0);

extern void warning(const char *msg, ...)
	       		__cold __format(printf, 1, 2);

extern double rnorm(double mean, double sd);

extern void *new_matrix(size_t dimc, ...) __cold;

extern void delete_matrix(void *matrix, size_t dimc, ...) __cold;

extern void make_matrix_row_cumulative(double row[], size_t len);

extern void make_2D_matrix_cumulative(double **matrix, size_t width,
				      size_t height);

extern void make_3D_matrix_cumulative(double ***matrix, size_t width,
				      size_t height, size_t depth);

extern bool prepare_uint64_probability_array(uint64_t array[], size_t len,
					     size_t extra_idx);


extern const char *timestamp;

extern void set_timestamp() __cold;
static inline const char *get_timestamp()
{
	return timestamp;
}

extern const char *command_line;

extern void set_command_line(int argc, const char * const * argv) __cold;
static inline const char *get_command_line()
{
	return command_line;
}

char *xstrdup(const char *s);
char *xstrdup2(const char *s1, const char *s2);
void *xmalloc(size_t size);
void *xrealloc(void *ptr, size_t size);

enum OutputType {
	DEFAULT = 0,
	TEXT,
	GZIP,
};


template <typename T>
static inline T max_array(const T array[], size_t len)
{
	T max = 0.0;
	for (size_t i = 0; i < len; i++)
		if (array[i] > max)
			max = array[i];
	return max;
}

template <typename T>
static inline T sum_array(const T array[], size_t len)
{
	T sum = 0;
	for (size_t i = 0; i < len; i++)
		sum += array[i];
	return sum;
}

template <typename T>
static inline void zero_array(T array[], size_t len)
{
	for (size_t i = 0; i < len; i++)
		array[i] = (T) 0;
}

static inline bool in_unit_interval(double n)
{
	return (n >= 0.0 && n <= 1.0);
}

static inline const char *bool_to_str(bool b)
{
	if (b)
		return "yes";
	else
		return "no";
}

#define ARRAY_LEN(array)      (sizeof(array) / sizeof((array)[0]))
#define ZERO_ARRAY(array)     (memset((array), 0, sizeof(array)))

extern const char dna_char_to_bin_tab[];

static inline unsigned char dna_char_to_bin(char c)
{
	return dna_char_to_bin_tab[((unsigned char)c)];
}

extern const char dna_bin_to_char_tab[];
static inline char dna_bin_to_char(unsigned char base)
{
	return dna_bin_to_char_tab[base];
}

static inline bool is_valid_dna_char(char c)
{
	return dna_char_to_bin(c) < 4;
}

static inline bool seq_contains_non_dna_chars(const char *seq, size_t len)
{
	while (len--)
		if (!is_valid_dna_char(*seq++))
			return true;
	return false;
}

extern const char dna_char_complement_tab[];
static inline char dna_char_complement(char c)
{
	return dna_char_complement_tab[(unsigned char)c];
}

/* Reverse complement an ASCII DNA sequence. */
extern void reverse_complement(char *p, size_t len);

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
extern size_t trim(char *s, size_t len);

#endif
