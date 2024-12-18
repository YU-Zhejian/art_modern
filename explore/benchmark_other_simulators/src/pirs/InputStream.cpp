#include "InputStream.h"
#include "util.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* See InputStream.h for documentation. */

enum OutputType InputStream::file_type_probe(const char *filename)
{
	const char *p = strrchr(filename, '.');
	if (p && strcmp(p, ".gz") == 0)
		return GZIP;
	else
		return TEXT;
}

int InputStream::gzip_status()
{
	const char *error_str;
	int errnum;
	if (gzeof(gzip.s_fp))
		return -1;
	error_str = gzerror(gzip.s_fp, &errnum);
	if (errnum == Z_ERRNO)
		fatal_error_with_errno("Error reading \"%s\"", s_filename);
	else
		fatal_error("zlib error while reading \"%s\": %s",
			    s_filename, error_str);
}

int InputStream::text_status()
{
	if (ferror(text.s_fp))
		fatal_error_with_errno("Error reading \"%s\"", s_filename);
	return -1;
}

void InputStream::open(const char *filename)
{
	bool is_stdin;
	if (strcmp(filename, "-") == 0) {
		is_stdin = true;
		s_filename = xstrdup("(stdin)");
		s_type = TEXT;
	} else {
		is_stdin = false;
		s_type = file_type_probe(filename);
		s_filename = xstrdup(filename);
	}
	//info("Opening input file \"%s\"\n", s_filename);

	switch (s_type) {
	case TEXT:
		if (is_stdin)
			text.s_fp = stdin;
		else
			text.s_fp = fopen(s_filename, "rb");
		if (!text.s_fp)
			goto fail;
		break;
	case GZIP:
		gzip.s_fp = gzopen(s_filename, "rb");
		if (!gzip.s_fp)
			goto fail;
		break;
	default:
		assert(0);
	}
	return;
fail:
	fatal_error_with_errno("Failed to open the file \"%s\"",
			       s_filename);
}

ssize_t InputStream::getline(char **lineptr, size_t *n)
{
	ssize_t ret;
	switch (s_type) {
	case TEXT:
		ret = ::getline(lineptr, n, text.s_fp);
		if (ret == -1)
			return text_status();
		break;
	case GZIP: {
			size_t offset = 0;
			if (!*lineptr) {
				*n = 128;
				*lineptr = (char*)xmalloc(*n);
			}
			while (1) {
				char *line = *lineptr;
				if (gzgets(gzip.s_fp, line + offset, *n - offset)) {
					ret = strlen(line);
					if (line[ret - 1] == '\n')
						return ret;
				} else {
					return gzip_status();
				}
				offset = *n - 1;
				*n *= 2;
				*lineptr = (char*)xrealloc(line, *n);
			}
			break;
		}
	default:
		assert(0);
	}
	return ret;
}

int InputStream::getc()
{
	int c;
	switch (s_type) {
	case TEXT:
		c = fgetc(text.s_fp);
		if (c == EOF)
			return text_status();
		break;
	case GZIP:
		c = gzgetc(gzip.s_fp);
		if (c == EOF)
			return gzip_status();
		break;
	default:
		assert(0);
	}
	return c;
}

void InputStream::ungetc(char c)
{
	switch (s_type) {
	case TEXT:
		::ungetc(c, text.s_fp);
		break;
	case GZIP:
		gzungetc(c, gzip.s_fp);
		break;
	default:
		assert(0);
	}
}

void InputStream::close()
{
	if (!s_filename)
		return;
	//info("Closing input file \"%s\"\n", s_filename);
	switch (s_type) {
	case TEXT:
		if (text.s_fp)
			fclose(text.s_fp);
		break;
	case GZIP:
		if (gzip.s_fp)
			gzclose(gzip.s_fp);
		break;
	default:
		assert(0);
	}
	free(s_filename);
	s_filename = NULL;
}

InputStream::~InputStream()
{
	if (s_type != DEFAULT)
		this->close();
}
